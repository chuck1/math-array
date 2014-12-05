CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(CMAKE_INSTALL_PREFIX /scratch/usr)

PROJECT(math-array)

FIND_PACKAGE(CMakeHelper)
FIND_PACKAGE(galaxy_log)



cmh_header_only()

ADD_SUBDIRECTORY(test)

