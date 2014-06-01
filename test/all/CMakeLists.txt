project(test-all)

FIND_PACKAGE(CMakeHelper)
FIND_PACKAGE(math-array)


SET(BOOST_ROOT /scratch/usr)
SET(BOOST_LIBRARY_DIR /scratch/usr/lib)
SET(Boost_NO_SYSTEM_PATHS ON)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
FIND_PACKAGE(
	Boost 1.55
	COMPONENTS log system thread
	REQUIRED)


SET(libs ${libs}
	Galaxy-Log
	${Boost_LIBRARIES}
	pthread
	rt
	)




link_exe()



