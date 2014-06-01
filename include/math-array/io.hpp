#ifndef MATH_ARRAY_IO_HPP
#define MATH_ARRAY_IO_HPP
/*
template<typename T> std::ostream&	operator<<(std::ostream& strm, std::vector<T> const & vec) {
	auto width = strm.width();
	for(T t : vec) {
		strm.width(width);
		strm << t;
	}
	return strm;
}
*/
/*
template<typename STREAM, typename T> STREAM&	operator<<(STREAM& strm, std::vector<T> const & vec) {
	auto width = strm.width();
	for(T t : vec) {
		strm.width(width);
		strm << t;
	}
	return strm;
}
*/
/*
template< typename CharT, typename TraitsT, typename T > inline
std::basic_ostream< CharT, TraitsT >& operator<< (std::basic_ostream< CharT, TraitsT >& strm, std::vector<T> const & vec) {
	auto width = strm.width();
	for(T t : vec) {
		strm.width(width);
		strm << t;
	}
	return strm;
}
*/

template< typename CharT, typename TraitsT > inline
std::basic_ostream< CharT, TraitsT >& operator<< (std::basic_ostream< CharT, TraitsT >& strm, std::vector<size_t> const & vec) {
	auto width = strm.width();
	for(auto t : vec) {
		strm.width(width);
		strm << t;
	}
	return strm;
}

#endif
