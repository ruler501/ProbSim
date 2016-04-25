#import <algorithm>
#import <functional>
#import <iomanip>
#import <iostream>
#import <map>
#import <numeric>
#import <sstream>
#import <string>
#import <thread>
#import <type_traits>
#import <utility>
#import <vector>

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
            if (it->second <= 0) dist.erase(it++);
            else it++;
        }
    }

    ProbabilityDist(const std::map<T, double>& d)
    : dist(d) {prune();}
    
    template <typename U>
    friend class recDefD2D;
    
    template <typename U>
    friend class recDefTD;

    template <typename Callable, typename U, typename... Args>
    friend auto distToDists(Callable tr, ProbabilityDist<U> p, Args... args) -> typename std::result_of<Callable(U, typename extractT<Args>::type...)>::type; 
    
    template <typename Callable, typename U, typename... Args>
    friend auto transformDists(Callable tr, ProbabilityDist<U> p, Args... args) -> ProbabilityDist<typename std::result_of<Callable(U, typename extractT<Args>::type...)>::type>; 
    
    friend ProbabilityDist<int> sumProb(int n, int s);

    template <typename U>
    friend class ProbabilityDist;
public:
    ProbabilityDist()
    : dist{{T(),1}} {}

    ProbabilityDist(T start)
    : dist{{start,1}} {}

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

ProbabilityDist<int> sumProb(int n, int s){
    if(n <= 0 || s <= 0) return ProbabilityDist<int>(0);
    double iProb = 1.f/s;
    std::map<int,double> distN;
    for(int i=1; i <= s; i++) distN[i] = iProb;

    for(int i=2; i <= n; i++){
        std::map<int, double> tDist;
        int mVal = s*i;
        for(int j=i; j <= mVal; j++){
            double tVal = 0;
            int imVal = mVal - s;
            for(int k=i-1; k <= imVal; k++)
                if(1 <= j-k && j - k <= s)
                    tVal += distN[k];
            tDist[j] = tVal/s;
        }
        distN = tDist;
    }
    return ProbabilityDist<int>(distN);
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
    int threads = std::thread::hardware_concurrency();
    std::vector<std::thread*>runThreads;
    std::vector<std::map<target, double>> dists(threads);
    std::vector<std::map<T, double>> sDists(threads);
    int count = 0;
    for(auto& e : p.dist){
       sDists[count++ % threads][e.first] = e.second;
    }
    for(int i = 0; i < threads; i++){
        if(!sDists[i].empty()){
            runThreads.push_back(new std::thread(&(recDefD2D<target>::template call<decltype(curry(tr)), T, Args...>), curry(tr), std::ref(dists[i]), std::vector<double>(), ProbabilityDist<T>(sDists[i]), args...));
        }
    }
    std::map<target,double> dist;
    for(int i=0; i < runThreads.size(); i++){
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
    int threads = std::thread::hardware_concurrency();
    std::vector<std::thread*>runThreads;
    std::vector<std::map<target, double>> dists(threads);
    std::vector<std::map<T, double>> sDists(threads);
    int count = 0;
    for(auto& e : p.dist){
       sDists[count++ % threads][e.first] = e.second;
    }
    for(int i = 0; i < threads; i++){
        if(!sDists[i].empty()){
            runThreads.push_back(new std::thread(&(recDefTD<target>::template call<decltype(curry(tr)), T, Args...>), curry(tr), std::ref(dists[i]), std::vector<double>(), ProbabilityDist<T>(sDists[i]), args...));
        }
    }
    std::map<target,double> dist;
    for(int i=0; i < runThreads.size(); i++){
            runThreads[i]->join();
            for(auto& e : dists[i]) dist[e.first] = getWithDef(dist, e.first, (double)0) + e.second;
    }

    return ProbabilityDist<target>(dist);
}

double sumAny(int a, float b, double c){
    return a + b + c;
}

ProbabilityDist<int> rollAgain(int n, int s){
    return sumProb(n-1,s);
}

bool heroWon(std::pair<int,int> life){
    return life.first > 0 && life.second <= 0;
}

bool battleDone(std::pair<int, int> life){
    return life.first <= 0 || life .second <= 0;
}

std::pair<int,int> newLife(std::pair<int, int> life, int hd, int bd){
    if (life.first > 0){
        life.second = std::max(0, life.second - hd);
    }
    if(life.second > 0){
        life.first = std::max(0, life.first - bd);;
    }
    return life;
}
    

int main(int argc, char* argv[]){
    //ProbabilityDist<int> d6 = sumProb(1,200);
    //ProbabilityDist<float> d6f(d6);
    //ProbabilityDist<double> d6d(d6);
    /* ProbabilityDist<int> d20 = sumProb(1,20); */
    /* //std::cout << d6 << std::endl << d6f << std::endl << d6d << std::endl; */
    //std::cout << ProbabilityDist<int>(transformDists(sumAny, d6, d6f, d6d)) << std::endl << std::endl;
    /* std::cout << sumProb(0,6) << std::endl << std::endl; */
    //std::cout << distToDists(rollAgain, sumProb(1,1000), sumProb(1,2)) << std::endl;
    ProbabilityDist<int> hd = sumProb(10,4);
    ProbabilityDist<int> bd = sumProb(2,20);
    ProbabilityDist<std::pair<int,int>> life(std::make_pair(100,100));
    ProbabilityDist<bool> won(false);
    ProbabilityDist<bool> done(false);
    while(done.valueAt(false) > 0.0001){
         life = transformDists(newLife, life, hd, bd);
         won = transformDists(heroWon, life);
         done = transformDists(battleDone, life);
         std::cout << won << '\t' << done << '\t';
         std::cout.unsetf(std::ios_base::fixed);
         std::cout << done.valueAt(false) << std::endl;
    }
}


//template <typename T, typename... Args>
//void plotDists(ProbabilityDist<T>, Args... args);
