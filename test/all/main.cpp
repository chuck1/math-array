#include <math-array/array.hpp>

int main(int ac, char** av) {
	
	auto arr1 = make_uninit<int,2>({3,4});
	//auto arr2 = make_zeros<int,3>({1,2,3});
	
	for(size_t i : range(2)) {
		for(size_t j : range(3)) {
			//for(size_t k : range(4)) {
				arr1->get(i,j) = i + j;
				std::cout << arr1->get(i,j) << std::endl;
			//}
		}
	}
	
	auto arr2 = arr1->sub({0,0},{2,3});
}





