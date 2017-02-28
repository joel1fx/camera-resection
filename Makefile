
LIB_OBJS=function_object.o

CXX_FLAGS=-Wall -O2

camera_resection: camera_resection.o $(LIB_OBJS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

function_object.o: function_object.cpp function_object.h
	$(CXX) $(CXX_FLAGS) -c function_object.cpp

camera_resection.o: camera_resection.cpp function_object.h
	$(CXX) $(CXX_FLAGS) -c camera_resection.cpp

clean:
	rm *.o
