#ifndef _PDIST_H_
#define _PDIST_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

//#include "Cache.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/*****************************************************************************
 * BEGIN BOILERPLATE UTILITIES                                               *
 *****************************************************************************/
template <std::size_t... Indices>
struct indices {
    using next = indices<Indices..., sizeof...(Indices)>;
};
template <std::size_t N>
struct build_indices {
    using type = typename build_indices<N-1>::type::next;
};
template <>
struct build_indices<0> {
    using type = indices<>;
};
template <std::size_t N>
using BuildIndices = typename build_indices<N>::type;

template <size_t num_args>
struct unpack
{
private:
    template <typename FuncType, size_t... I>
    auto call(FuncType &f, std::vector<int> &args, indices<I...>)
      -> decltype( f(args[I]...) )
    {
        return f(args[I]...);
    }

public:
    template <typename FuncType>
    auto operator () (FuncType &f, std::vector<int> &args)
      -> decltype( std::declval<unpack<num_args>>().call(f, args, BuildIndices<num_args>{}) )
    {
        return call(f, args, BuildIndices<num_args>{});
    }
}; 

namespace _dtl {

    template <typename FUNCTION> struct
    _curry;

    // specialization for functions with a single argument
    template <typename R,typename T> struct
    _curry<std::function<R(T)>> {
        using
        type = std::function<R(T)>;

        const type
        result;

        _curry(type fun) : result(fun) {}

    };

    // recursive specialization for functions with more arguments
    template <typename R,typename T,typename...Ts> struct
    _curry<std::function<R(T,Ts...)>> {
        using
        remaining_type = typename _curry<std::function<R(Ts...)> >::type;

        using
        type = std::function<remaining_type(T)>;

        const type
        result;

        _curry(std::function<R(T,Ts...)> fun)
        : result (
            [=](const T& t) {
                return _curry<std::function<R(Ts...)>>(
                    [=](const Ts&...ts){ 
                        return fun(t, ts...); 
                    }
                ).result;
            }
        ) {}
    };
}

template <typename R,typename...Ts> auto
curry(const std::function<R(Ts...)> fun)
-> typename _dtl::_curry<std::function<R(Ts...)>>::type
{
    return _dtl::_curry<std::function<R(Ts...)>>(fun).result;
}

template <typename R,typename...Ts> auto
curry(R(* const fun)(Ts...))
-> typename _dtl::_curry<std::function<R(Ts...)>>::type
{
    return _dtl::_curry<std::function<R(Ts...)>>(fun).result;
}

template <template<class,class,class...> class C, typename K, typename V, typename... Args>
V getWithDef(const C<K,V,Args...>& m, K const& key, const V & defval)
{
    typename C<K,V,Args...>::const_iterator it = m.find( key );
    if (it == m.end())
        return defval;
    return it->second;
}
/*****************************************************************************
 * END BOILERPLATE UTILITIES                                                 *
 *****************************************************************************/

template <typename T>
struct extractT{};

template <typename T>
class ProbabilityDist;

template<typename T>
struct extractT<ProbabilityDist<T>>{
    typedef T type;
};

template <typename T>
class ProbabilityDist{
protected:
    std::map<T,double> dist;

    void prune(){
        return;
        for(auto it=dist.begin(); it != dist.end();){
            if (it->second <= 0) it = dist.erase(it);
            else it++;
        }
    }

    template <typename U>
    friend struct recDefD2D;
    
    template <typename U>
    friend struct recDefTD;

    template <typename Callable, typename U, typename... Args>
    friend auto distToDists(Callable tr, ProbabilityDist<U> p, Args... args) -> typename std::result_of<Callable(U, typename extractT<Args>::type...)>::type; 
    
    template <typename Callable, typename U, typename... Args>
    friend auto transformDists(Callable tr, ProbabilityDist<U> p, Args... args) -> ProbabilityDist<typename std::result_of<Callable(U, typename extractT<Args>::type...)>::type>; 
    
    friend ProbabilityDist<int> sumProb(int n, int s);

    friend void plotDist(ProbabilityDist<T> d, const char* format=nullptr){
        std::vector<T> x(d.dist.size());
        std::vector<double> y(d.dist.size());
        int i = 0;
        for (auto& e : d.dist){
            x[i] = e.first;
            y[i] = e.second;
            i++;
        }
        if(format != nullptr) plt::plot(x, y, format);
        else plt::plot(x,y);
    }

    template <typename U>
    friend class ProbabilityDist;
public:
    ProbabilityDist()
    : dist{{T(),1}} {}

    ProbabilityDist(T start)
    : dist{{start,1}} {}

    ProbabilityDist(const std::map<T, double>& d)
    : dist(d) {prune();}
    
    std::map<T,double>& getDist(){ return dist; }

    template <typename U>
    ProbabilityDist(ProbabilityDist<U> d){
        for(auto e : d.dist){
            dist[U(e.first)] = e.second;
        }
    }
//    T average();

//    T mVal();

//    std::vector<T> mVal();

    friend std::ostream& operator<<(std::ostream& os, const ProbabilityDist<T> dist){
        auto it = dist.dist.begin();
        os << std::fixed << std::setprecision(2);
        os << "[(" << it->first << ", " << it->second*100 << "%)";
        it++;
        for(; it != dist.dist.end(); it++){
            os << "; (" << it->first << ", " << it->second*100 << "%)";
        }
        return os << "]";
    }

//    ProbabilityDist<T> cumulativeDistribution();

    double valueAt(T val){
        auto iter = dist.find(val);
        if(iter == dist.end()) return 0;
        else return iter->second;
    }
};

std::map<int, double> addDist(std::map<int, double> x,
                              std::map<int, double> y){
    std::map<int, double> tDist;
    int minX = x.begin()->first;
    int maxX = x.rbegin()->first;
    int minY = y.begin()->first;
    int maxY = y.rbegin()->first;
    int maxVal = maxX + maxY;
    int minVal = minX + minY;
    for(int i=minVal; i <= maxVal; i++){
        double tVal = 0;
        for(int j=minX; j <= maxX; j++){
            tVal += x[j]*y[i-j];
        }
        tDist[i] = tVal;
    }
    return tDist;
}

ProbabilityDist<int> sumProb(int n, int s){
    std::map<int, double> distX;
    std::map<int, double> distZ;
    distZ[0] = 1.f;
    double baseProb = 1.f/s;
    for(int i=1; i<=s; i++) { distX[i] = baseProb; }
    while(n > 1){
        if(n % 2){
            distZ = addDist(distX, distZ);
            distX = addDist(distX, distX);
        }
        else{
            distX = addDist(distX, distX);
        }
        n /= 2;
    }
    return ProbabilityDist<int>(addDist(distX, distZ));
}

template <typename T>
struct recDefD2D{
    static void call(ProbabilityDist<T> res, std::map<T, double>& dist, std::vector<double>& pargs){
        double val = std::accumulate(pargs.begin(), pargs.end(), (double)1, std::multiplies<double>());
        for(auto& e: res.dist){
            if(e.second > 0)
                dist[e.first] = getWithDef(dist, e.first, (double)0) + val*e.second;
        }
    }

    template<typename Callable, typename U, typename... RArgs>
    static void call(Callable tr, std::map<T, double>& dist, std::vector<double> pargs, ProbabilityDist<U> carg, RArgs... rargs){
        for(auto& e : carg.dist){
            if(e.second > 0){
                std::vector<double> npargs(pargs);
                npargs.push_back(e.second);
                recDefD2D<T>::call(tr(e.first), dist, npargs, rargs...);
            }
        }
    }
};

template <typename Callable, typename T, typename... Args>
auto distToDists(Callable tr, ProbabilityDist<T> p, Args... args) -> typename std::result_of<Callable(T, typename extractT<Args>::type...)>::type{ 
    using target = typename extractT<typename std::result_of<Callable(T, typename extractT<Args>::type...)>::type>::type;
    using ftype = std::function<ProbabilityDist<target>(T, typename extractT<Args>::type...)>;
    int threads = std::thread::hardware_concurrency();
#ifdef DISABLEMULTITHREAD
    threads = 1;
#endif//DISABLEMULTITHREAD
    std::vector<std::thread*>runThreads;
    std::vector<std::map<target, double>> dists(threads);
    std::vector<std::map<T, double>> sDists(threads);
    int count = 0;
    for(auto& e : p.dist){
       sDists[count++ % threads][e.first] = e.second;
    }
    for(int i = 0; i < threads; i++){
        if(!sDists[i].empty()){
            runThreads.push_back(new std::thread(
                        recDefD2D<target>::template call<decltype(curry(ftype(tr))), T, Args...>,
                        curry(ftype(tr)),
                        std::ref(dists[i]),
                        std::vector<double>(),
                        ProbabilityDist<T>(sDists[i]),
                        args...));
        }
    }
    std::map<target,double> dist;
    for(unsigned int i=0; i < runThreads.size(); i++){
            runThreads[i]->join();
            for(auto& e : dists[i]) dist[e.first] = getWithDef(dist, e.first, (double)0) + e.second;
    }

    return ProbabilityDist<target>(dist);
}

template <typename T>
struct recDefTD{
    static void call(T res, std::map<T, double>& dist, std::vector<double>& pargs){
        double val = std::accumulate(pargs.begin(), pargs.end(), (double)1, std::multiplies<double>());
        dist[res] = getWithDef(dist, res, (double)0) + val;
    }

    template<typename Callable, typename U, typename... RArgs>
    static void call(Callable tr, std::map<T, double>& dist, std::vector<double> pargs, ProbabilityDist<U> carg, RArgs... rargs){
        for(auto& e : carg.dist){
            if(e.second > 0){
                std::vector<double> npargs(pargs);
                npargs.push_back(e.second);
                recDefTD<T>::call(tr(e.first), dist, npargs, rargs...);
            }
        }
    }
};

template <typename Callable, typename T, typename... Args>
auto transformDists(Callable tr, ProbabilityDist<T> p, Args... args) -> ProbabilityDist<typename std::result_of<Callable(T, typename extractT<Args>::type...)>::type>{ 
    using target = typename std::result_of<Callable(T, typename extractT<Args>::type...)>::type;
    using ftype = std::function<target(T, typename extractT<Args>::type...)>;
    int threads = std::thread::hardware_concurrency();
#ifdef DISABLEMULTITHREAD
    threads = 1;
#endif//DISABLEMULTITHREAD
    std::vector<std::thread*>runThreads;
    std::vector<std::map<target, double>> dists(threads);
    std::vector<std::map<T, double>> sDists(threads);
    int count = 0;
    for(auto& e : p.dist){
       sDists[count++ % threads][e.first] = e.second;
    }
    for(int i = 0; i < threads; i++){
        if(!sDists[i].empty()){
            runThreads.push_back(new std::thread(
                        recDefD2D<target>::template call<decltype(curry(ftype(tr))), T, Args...>,
                        curry(ftype(tr)),
                        std::ref(dists[i]),
                        std::vector<double>(),
                        ProbabilityDist<T>(sDists[i]),
                        args...));
        }
    }
    std::map<target,double> dist;
    for(unsigned int i=0; i < runThreads.size(); i++){
            runThreads[i]->join();
            for(auto& e : dists[i]) dist[e.first] = getWithDef(dist, e.first, (double)0) + e.second;
    }

    return ProbabilityDist<target>(dist);
}

template<typename Numeric>
double averageDist(ProbabilityDist<Numeric> d){
    double t = 0;
    for(auto& e : d.getDist()){
        t += e.first*e.second;
    }
    return t;
}

#endif
