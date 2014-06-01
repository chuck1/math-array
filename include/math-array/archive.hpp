#ifndef MATH_ARRAY_ARCHIVE_HPP
#define MATH_ARRAY_ARCHIVE_HPP

#include <iostream>
#include <fstream>
#include <vector>

namespace math {
	class basic_binary_oarchive {
		public:
			basic_binary_oarchive(std::ofstream& ofs): ofs_(ofs) {
			}
	
			
			template<typename T> basic_binary_oarchive&		operator<<(T&& t) {
				T c(t);
				ofs_.write((char*)&c, sizeof(T));
				return *this;
			}
			template<typename T> basic_binary_oarchive&		operator<<(T const & t) {
				ofs_.write((char*)&t, sizeof(T));
				return *this;
			}
			template<typename T> basic_binary_oarchive&		operator<<(std::vector<T>& t) {
				operator<<(t.size());
				for(auto it = t.cbegin(); it != t.cend(); ++it) {
					operator<<(*it);
				}
				return *this;
			}

			std::ofstream&		ofs_;
	};
}



#endif



