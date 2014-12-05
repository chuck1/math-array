#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>

#include <gal/log/log.hpp>

#include <math-array/io.hpp>
#include <math-array/range.hpp>
#include <math-array/archive.hpp>


/** @file array.hpp
*/

typedef math::basic_binary_oarchive oarchive;


namespace math {


	template<typename T, int N> class __array;

	/** @typedef array
	*/
	template<typename T, int N> using array = std::shared_ptr< __array<T,N> >;

#define FOR(i,level) for(size_t i = 0; i < n_[level]; i++)

	template <int D, typename U> struct __multivec {
		typedef std::vector< typename __multivec<D-1,U>::vec_type > vec_type;
	};
	template <typename U> struct __multivec<1,U> {
		typedef std::vector<U> vec_type;
	};

	template<int D, typename U> using multivec = typename __multivec<D,U>::vec_type;


	template <int D, typename U> struct Initializer_list {
		typedef std::initializer_list<typename Initializer_list<D-1,U>::list_type > list_type;
	};
	template <typename U> struct Initializer_list<1,U> {
		typedef std::initializer_list<U> list_type;
	};



	template<typename T, int N> array<T,N>		make_array();
	template<typename T, int N> array<T,N>		make_array(std::vector<size_t>, std::initializer_list<T> il);
	template<typename T, int N> array<T,N>		make_ones(array<size_t,1> v);
	template<typename T, int N> array<T,N>		make_ones(std::vector<size_t> v);
	template<typename T, int N> array<T,N>		make_uninit(std::vector<size_t> n);






	template <int N, typename T> std::vector<int>	size(typename Initializer_list<N,T>::list_type const & il) {
		std::vector<int> s;

		s.push_back(il.size());

		size(s, il[0]);

		return s;
	}


	template<typename T> void print(T* t, int n) {
		for(int i = 0; i < n; ++i) {
			std::cout << std::setw(16) << t[i];
		}
		std::cout << std::endl;
	}
	template<typename T> void print(std::vector<T> v) {
		for(T i : v) {
			std::cout << std::setw(16) << i;
		}
		std::cout << std::endl;
	}


	namespace arr {
		template<typename...> int prod(int p) {
			return p;
		}
		template<typename... B> int prod(int p, int a, B... b) {
			return prod(p*a,b...);
		}
	}


	template<typename...> void fill(int*) {
	}
	template<int A, int... B> void fill(int* i) {
		*i = A;
		fill<B...>(i+1);
	}


	template<typename...> void cum(int*, int) {
	}
	template<int A, int... B> void cum(int* i, int s) {
		*i = s/A;
		cum<B...>(i+1,*i);
	}

	template<typename T, int N> class __array: public std::enable_shared_from_this< __array<T,N> > {
		public:
			typedef std::shared_ptr< __array<T,N> >			shared;
			typedef typename Initializer_list<N,T>::list_type	init_list;
		public:
			__array():
				size_(0),
				v_(0)
		{
		}
			__array(__array<T,N> const & rhs):
				size_(0),
				v_(0)
		{
			alloc(rhs.n_);
			memcpy(v_, rhs.v_, size_ * sizeof(T));
		}
			/** @brief move constructor
			*/
			__array(__array<T,N> && rhs) {
				std::swap(c_, rhs.c_);
				std::swap(n_, rhs.n_);
				std::swap(size_, rhs.size_);
				std::swap(v_, rhs.v_);
			}
			/** @brief move assignment
			*/
			__array<T,N>&		operator=(__array<T,N> && rhs) {
				std::swap(c_, rhs.c_);
				std::swap(n_, rhs.n_);
				std::swap(size_, rhs.size_);
				std::swap(v_, rhs.v_);
				return *this;
			}
		public:
			/** @name Allocation
			 * @{
			 */
			void					alloc(math::array<int,1> shape_arr) {
				std::vector<size_t> shape;
				for(int i : *shape_arr) shape.push_back(i);
				alloc(shape);
			}
			void					alloc(math::array<size_t,1> shape_arr) {
				std::vector<size_t> shape;
				for(int i : *shape_arr) shape.push_back(i);
				alloc(shape);
			}
			void					alloc(std::vector<size_t> n) {
				assert(this);
				assert(n.size() == N);

				n_ = n;

				size_ = 1;
				for(size_t i : n_) {
					assert(i > 0);
					size_ *= i;
				}

				c_.clear();
				size_t s(size_);
				for(size_t i : n_) {
					assert(i > 0);
					//std::cout << "i s (i > 0)" << i << " " << s << " " << (i>0) << std::endl;
					s /= i;
					c_.push_back(s);
				}

				v_ = new T[size_];

				//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", debug) << "n    " << n_ << std::endl;
				//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", debug) << "c    " << c_ << std::endl;
				//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", debug) << "size " << size_ << std::endl;
			}
			/** @} */
			void					zeros() {
				for(size_t i = 0; i < size_; ++i) {
					v_[i] = 0;
				}
			}
			void					ones() {
				for(size_t i = 0; i < size_; ++i) {
					v_[i] = 1;
				}
			}

			void					set(std::initializer_list<T> src) {
				__set(v_, src);
			}
			void					set(std::vector<int> i, std::initializer_list<T> src) {
				T* dst = __get(v_, i, c_);
				__set(dst, src);
			}
			void					__set(T* dst, std::initializer_list<T> src) {
				for(T s : src) {
					*dst = s;
					dst++;
				}
			}

			void					set(std::vector<int> i, T* src, int len) {
				T* dst = __get(v_, i, c_);

				memcpy(dst, src, len * sizeof(T));
			}


			shared					transpose() {
				if(N != 2) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "Transpose for array with N==" << N << " not implemented" << std::endl;
					abort();
				}

				auto ret = make_uninit<T,N>({n_[1], n_[0]});

				FOR(i, 0) {
					FOR(j, 1) {
						ret->get(j,i) = get(i,j);
					}
				}

				return ret;
			}
			shared					transpose_self() {
				if(N != 2) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "Transpose for array with N==" << N << " not implemented" << std::endl;
					abort();
				}

				auto ret = transpose();

				std::swap(*this,*ret);

				return __array<T,N>::shared_from_this();
			}
			shared					rot90(int i) {
				if(N != 2) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "rot90 for N==" << N << " not implemented" << std::endl;
					abort();
				}
				if(i != 1) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "rot90 for i==" << i << " not implemented" << std::endl;
					abort();
				}

				auto ret = make_uninit<T,N>({n_[1], n_[0]});

				FOR(i, 0) {
					FOR(j, 1) {
						ret->get(n_[1] - 1 - j,i) = get(i,j);
					}
				}

				return ret;
			}
			shared					rot90_self(int i) {
				if(N != 2) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "rot90 for N==" << N << " not implemented" << std::endl;
					abort();
				}
				if(i != 1) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "rot90 for i==" << i << " not implemented" << std::endl;
					abort();
				}

				auto ret = rot90(i);

				std::swap(*this,*ret);

				return __array<T,N>::shared_from_this();
			}
			shared					fliplr() {
				if(N != 2) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "fliplr for N==" << N << " not implemented" << std::endl;
					abort();
				}

				auto ret = std::make_shared< __array<T,N> >(*this);

				ret->fliplr_self();

				return ret;
			}
			shared					fliplr_self() {
				if(N != 2) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "fliplr for N==" << N << " not implemented" << std::endl;
					abort();
				}
				size_t half = (n_[1] - (n_[1] % 2)) / 2;
				FOR(i, 0) {
					for(size_t j = 0; j < half; ++j) {
						std::swap(get(i,j), get(i, n_[1] - 1 - j));
					}
				}

				return __array<T,N>::shared_from_this();
			}
			std::vector<T>				ravel() {
				std::vector<T> ret;
				for(T* p = v_; p < (v_ + size_); ++p) {
					ret.push_back(*p);
				}
				return ret;
			}
			/** @name access
			 * @{ */
			T&						operator()(int i) {
				i %= size_;
				return v_[i];
			}
			T&						get(int i[N]) {
				std::vector<int> iv(i, i+N-1);
				return *(__get(v_, c_, i, i+N-1));
			}
			T*						__get(T* t, int* c, int* i, int* i_last) {
				t += (*c) * (*i);

				if(i == i_last) return t;

				return __get(t,c+1,i+1,i_last);
			}
			/** @brief get
			 * std::vector
			 */
			T&						get(std::vector<int> const & i) {
				assert(i.size() <= N);
				std::vector<size_t> iu;
				for(size_t a = 0; a < N; ++a) {
					int b = (i[a] + n_[a]) % n_[a];
					assert((b >= 0) && ((size_t)b < n_[a]));
					iu.push_back((size_t)b);
				}
				return *(__get1(v_, c_.cbegin(), n_.cbegin(), iu.cbegin(), iu));
			}
			/** @brief get
			 * std::vector
			 */
			T&						get(std::vector<size_t> const & i) {
				assert(i.size() <= N);
				for(size_t a = 0; a < N; ++a) {
					assert(i[a] < n_[a]);
				}
				return *(__get1(v_, c_.cbegin(), n_.cbegin(), i.cbegin(), i));
			}
			T*						__get1(
					T* t,
					std::vector<size_t>::const_iterator c,
					std::vector<size_t>::const_iterator n,
					std::vector<size_t>::const_iterator it,
					std::vector<size_t> const & i) {
				if(it == i.cend()) return t;

				t += (*c) * (*it);

				return __get1(
						t,
						std::next(c),
						std::next(n),
						std::next(it),
						i);
			}
			/** @brief get
			 * variadic input
			 */
			template<typename... I> T&			get(I... b) {
				T* ptr = __get<I...>(v_, c_.begin(), n_.begin(), b...);

				//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", debug) << "get " << (int)(ptr - v_) << std::endl;

				return *(ptr);
			}
			template<typename... I> T*			__get(T* t, std::vector<size_t>::const_iterator c, std::vector<size_t>::const_iterator n, I... b) {
				return t;
			}
			template<typename A, typename... B> T*		__get(T* t, std::vector<size_t>::const_iterator c, std::vector<size_t>::const_iterator n, A a, B... b) {
				assert(c != c_.end());
				assert(n != n_.end());

				a = (a + (*n)) % (*n);
				assert((a >= 0) && ((size_t)a < (*n)));

				t += (*c) * a;

				c++;
				n++;

				return __get(t,c,n,b...);
			}
			void						copy(
					T* src,
					T* dst,
					T* const & src0,
					T* const & dst0,
					std::vector<size_t>::const_iterator i_src_b,
					std::vector<size_t>::const_iterator i_src_e,
					std::vector<size_t>::const_iterator i_dst_b,
					std::vector<size_t>::const_iterator i_dst_e,
					std::vector<size_t>::const_iterator c_src,
					std::vector<size_t>::const_iterator c_dst,
					std::vector<size_t> const & vec_i_src_b) {

				//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", debug)
				//	<< "copy *(src + " << (src-src0) << ") to *(dst + " << (dst-dst0) << ")" << std::endl;


				if(i_src_b == vec_i_src_b.cend()) {
					*dst = *src;
					return;
				}


				size_t i_src = *i_src_b;
				size_t i_dst = *i_dst_b;

				for(; i_src < *i_src_e; ++i_src, ++i_dst) {

					copy(
							src + i_src * (*c_src),
							dst + i_dst * (*c_dst),
							src0,
							dst0,
							std::next(i_src_b),
							std::next(i_src_e),
							std::next(i_dst_b),
							std::next(i_dst_e),
							std::next(c_src),
							std::next(c_dst),
							vec_i_src_b);


				}
			}



			/** @} */
			bool		equal_size(__array<T,N> const & rhs) {
				for(int i = 0; i < N; ++i) {
					if(n_[i] != rhs.n_[i]) return false;
				}
				return true;
			}

			__array<T,N>				operator+(__array<T,N> const & rhs) {
				__array<T,N> ret(*this);
				ret += rhs;
				return ret;
			}
			/** @name self arithmetic
			 * @{ */
			__array<T,N>&				operator+=(__array<T,N> const & rhs) {
				assert(equal_size(rhs));

				for(int i = 0; i < size_; ++i) {
					v_[i] += rhs.v_[i];
				}

				return *this;
			}
			__array<T,N>&				operator*=(__array<T,N> const & rhs) {
				assert(equal_size(rhs));

				for(int i = 0; i < size_; ++i) {
					v_[i] *= rhs.v_[i];
				}

				return *this;
			}
			__array<T,N>&				operator/=(__array<T,N> const & rhs) {
				assert(equal_size(rhs));

				for(int i = 0; i < size_; ++i) {
					v_[i] /= rhs.v_[i];
					assert(!std::isinf(v_[i]));
				}

				return *this;
			}
			/** @} */
			/** @name scalar self arithmatic @{ */
			__array<T,N>&				operator+=(T const & rhs) {
				for(int i = 0; i < size_; ++i) {
					v_[i] += rhs;
				}

				return *this;
			}
			__array<T,N>&				operator-=(T const & rhs) {
				for(int i = 0; i < size_; ++i) {
					v_[i] -= rhs;
				}

				return *this;
			}
			__array<T,N>&				operator*=(T const & rhs) {
				for(int i = 0; i < size_; ++i) {
					v_[i] *= rhs;
				}

				return *this;
			}
			__array<T,N>&				operator/=(T const & rhs) {
				for(int i = 0; i < size_; ++i) {
					v_[i] /= rhs;
				}

				return *this;
			}
			/** @} */
			/** @name shared self arithmetic
			 * @{ */
			shared					add_self(std::shared_ptr< __array<T,N> > const & rhs) {
				this->operator+=(*rhs);
				return __array<T,N>::shared_from_this();
			}
			shared					multiply_self(std::shared_ptr< __array<T,N> > const & rhs) {
				this->operator*=(*rhs);
				return __array<T,N>::shared_from_this();
			}
			shared					divide_self(std::shared_ptr< __array<T,N> > const & rhs) {
				this /= rhs;
				return *this;
			}
			/** @} */
			/** @name shared self scalar arithmetic
			 * @{ */
			shared					add_self(T const & rhs) {
				this->operator+=(rhs);
				return __array<T,N>::shared_from_this();
			}
			shared					multiply_self(T const & rhs) {
				operator*=(rhs);
				return __array<T,N>::shared_from_this();
			}
			shared					divide_self(T const & rhs) {
				operator/=(rhs);
				return __array<T,N>::shared_from_this();
			}
			/** @} */
			/** @name shared arithmetic
			 * @{ */
			shared					add(std::shared_ptr< __array<T,N> > const & rhs) {
				auto ret = std::make_shared< __array<T,N> >(*this);
				(*ret) += (*rhs);
				return ret;
			}
			shared					multiply(std::shared_ptr< __array<T,N> > const & rhs) {
				auto ret = std::make_shared< __array<T,N> >(*this);
				ret *= rhs;
				return ret;
			}
			shared					divide(std::shared_ptr< __array<T,N> > const & rhs) {
				auto ret = std::make_shared< __array<T,N> >(*this);
				ret /= rhs;
				return ret;
			}
			/** @} */
			/** @name shared scalar arithmetic
			 * @{ */
			shared					add(T const & rhs) {
				auto ret = std::make_shared< __array<T,N> >(*this);
				(*ret) += rhs;
				return ret;
			}
			shared					multiply(T const & rhs) {
				auto ret = std::make_shared< __array<T,N> >(*this);
				(*ret) *= rhs;
				return ret;
			}
			shared					divide(T const & rhs) {
				auto ret = std::make_shared< __array<T,N> >(*this);
				(*ret) /= rhs;
				return ret;
			}
			/** @} */
			/** @name math
			  @{ */
			shared								square() {

				auto ret = std::make_shared< __array<T,N> >(*this);
				ret->square_self();
				return ret;
			}
			template<typename U> std::shared_ptr< __array<U,N> >		ceil() {
				auto ret = make_uninit<U,N>(n_);
				for(size_t i = 0; i < size_; ++i) {
					ret->v_[i] = (U)(::ceil(v_[i]));
				}
				return ret;
			}
			shared								cumsum() {
				if(N != 1) {
					//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical) << "fliplr for N==" << N << " not implemented" << std::endl;
					abort();
				}

				auto ret = make_uninit<T,N>(n_);

				ret->get(0) = get(0);
				for(size_t i = 1; i < n_[0]; ++i) {
					ret->get(i) = ret->get(i-1) + get(i);
				}

				return ret;
			}
			void					square_self() {
				for(T* t = v_; t < (v_ + size_); ++t) {
					(*t) *= (*t);
				}
			}
			T					sum() {
				T s = 0;
				for(T* t = v_; t < (v_ + size_); ++t) {
					s += *t;
				}
				return s;
			}
			T					min() {
				T s = 1E37;
				for(T* t = v_; t < (v_ + size_); ++t) {
					if(*t < s) s = *t;
				}
				return s;
			}
			T					max() {
				T s = -1E37;
				for(T* t = v_; t < (v_ + size_); ++t) {
					if(*t > s) s = *t;
				}
				return s;
			}
			/** @name indexing @{ */
			std::vector<int>			offset2index(int o) {
				std::vector<int> ret;
				int a = o;
				for(int c : c_) {
					if(c == 1) {
						ret.push_back(a);
					} else {
						a = a % c;
						ret.push_back((o - a)/c);
					}
				}
				return ret;
			}
			/** @} */
			/** @name linear algebra
			 * @{ */
			std::shared_ptr< __array<T,N+1> >	gradient(array<T,N+1> d) {
				auto n = n_;
				n.push_back(N);

				auto ret = make_uninit<T,N+1>(n);

				//BOOST_LOG_CHANNEL_SEV(gal::log::lg, "math-array", critical)  << "sorry, not yet implemented" << std::endl;
				abort();

				//for(int i : range(size_)) {
				//}

				return ret;
			}
			T					fda_1_back(int i1, T h) {
				assert(N == 1);

				size_t iu0 = ((i1-1) + size_) % size_;
				size_t iu1 = (i1 + size_) % size_;

				return ((v_[iu1] - v_[iu0]) / h);
			}
			T					fda_2_back(int i2, T h) {
				assert(N == 1);

				size_t iu0 = ((i2-2) + size_) % size_;
				size_t iu1 = ((i2-1) + size_) % size_;
				size_t iu2 = (i2 + size_) % size_;

				return ((v_[iu2] - (2.0 * v_[iu1]) + v_[iu0])/(h*h));
			}
			/** @} */
			/** @name iterating
			 * @{ */
			T*					begin() {
				return v_;
			}
			T*					end() {
				return v_ + size_;
			}

			shared					sub(std::vector<int> beg_s, std::vector<int> end_s) {
				assert(beg_s.size() == end_s.size());

				std::vector<size_t> beg, end;

				// resolve negative indices
				for(size_t i = 0; i < end_s.size(); ++i) {
					beg.push_back((size_t)((beg_s[i] + n_[i]) % n_[i]));
					end.push_back((size_t)((end_s[i] + n_[i]) % n_[i]));
				}
				// calculate shape
				std::vector<size_t> shape = end;
				for(size_t i = 0; i < beg.size(); ++i) {
					shape[i] -= beg[i];
				}

				// alloc new array
				auto ret = std::make_shared< __array<T,N> >();
				ret->alloc(shape);

				std::vector<size_t> dst_b(N,0);

				// fill
				copy(
						v_,
						ret->v_,
						v_,
						ret->v_,
						beg.cbegin(),
						end.cbegin(),
						dst_b.cbegin(),
						ret->n_.cbegin(),
						c_.cbegin(),
						ret->c_.cbegin(),
						beg);
				return ret;
			}
			/** @} */
			/** @name info
			 * @{ */
			size_t					size() {
				return size_;
			}
			std::vector<size_t>			shape() {
				return n_;
			}
			/** @} */
			/** @name IO
			 * @{ */
			void					serialize(oarchive& ar, unsigned int const & version) {
				// assumes dimension and type are known
				ar << n_;
				for(size_t i = 0; i < size_; ++i) {
					ar << v_[i];
				}
			}
			/** @} */
		public://private:
			std::vector<size_t>		c_;
			std::vector<size_t>		n_;
			size_t				size_;
			T*				v_;
	};




	/*template<typename T, int N> array<T,N>		make_array() {
	  return std::make_shared< __array<T,N> >();
	  }*/

	template<typename T, int N> array<T,N>		make_copy(array<T,N> const & rhs) {
		auto arr = std::make_shared< __array<T,N> >(*rhs);
		return arr;
	}
	template<typename T, int N> array<T,N>		make_array_1(std::initializer_list<T> il) {
		std::vector<size_t> n({il.size()});
		return make_array<T,N>(n, il);
	}
	template<typename T, int N> array<T,N>		make_array(std::vector<size_t> n, std::initializer_list<T> il) {
		auto arr = std::make_shared< __array<T,N> >();
		arr->alloc(n);
		arr->set(il);
		return arr;
	}

	template<typename T, int N> math::array<T,N>		make_ones_arr(math::array<size_t,1> v) {
		auto arr = std::make_shared< __array<T,N> >();
		std::vector<size_t> a(std::begin(*v), std::end(*v));
		arr->alloc(a);
		arr->ones();
		return arr;
	}
	template<typename T, int N> math::array<T,N>		make_ones(std::vector<size_t> n) {
		auto arr = std::make_shared< __array<T,N> >();
		arr->alloc(n);
		arr->ones();
		return arr;
	}

	template<typename T, int N> math::array<T,N>		make_uninit(std::vector<size_t> n) {
		auto arr = std::make_shared< __array<T,N> >();
		arr->alloc(n);
		return arr;
	}

	template<typename T, int N> math::array<T,N>		make_zeros(std::vector<size_t> n) {
		auto arr = std::make_shared< __array<T,N> >();
		arr->alloc(n);
		arr->zeros();
		return arr;
	}

#define ARR_LOOP_START(arr) for(int p0 = 0; p0 < (arr->c_[0] * arr->n_[0]); p0 += arr->c_[0])

#define ARR_LOOP(arr, p, level) for(int p##level = p; p##level < (p + (arr->c_[level] * arr->n_[level])); p##level += arr->c_[level])


	template<typename T> array<T,1>		linspace(T s, T e, size_t n) {

		auto ret = make_uninit<T,1>({n});

		T step = (e - s) / ((T)(n - 1));

		for(int i : range(n)) {
			ret->get(i) = s + step * i;
		}

		return ret;
	}


	template<typename T> std::pair< array<T,2>, array<T,2> >		meshgrid(array<T,1> x, array<T,1> y) {

		auto nx = x->n_[0];
		auto ny = y->n_[0];

		auto X = make_uninit<T,2>({nx,ny});

		auto Y = make_uninit<T,2>({nx,ny});

		for(int i : range(nx)) {
			for(int j : range(ny)) {
				X->get(i,j) = x->get(i);
				Y->get(i,j) = y->get(j);
			}
		}

		return std::make_pair(X,Y);
	}

	}

#endif




