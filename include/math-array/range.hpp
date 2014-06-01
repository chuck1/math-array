#ifndef DIFF2D_MATH_HPP
#define DIFF2D_MATH_HPP

#include <vector>

template<typename T> std::vector<T>	range(T stop) {
	std::vector<T> ret;
	T step = 1;
	for(T i = 0; i < stop; i += step) {
		ret.push_back(i);
	}
	return ret;
}


#endif



