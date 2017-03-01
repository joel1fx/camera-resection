
#
# The MIT License (MIT)
#
# Copyright (c) 2015 Joel E. Merritt
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

LIB_OBJS=function_object.o g_n_solver.o

CXX_FLAGS=-Wall -O2

camera_resection: camera_resection.o $(LIB_OBJS)
	$(CXX) $(CXX_FLAGS) -o $@ $^

function_object.o: function_object.cpp function_object.h
	$(CXX) $(CXX_FLAGS) -c function_object.cpp

g_n_solver.o: g_n_solver.cpp g_n_solver.h
	$(CXX) $(CXX_FLAGS) -c g_n_solver.cpp

camera_resection.o: camera_resection.cpp function_object.h
	$(CXX) $(CXX_FLAGS) -c camera_resection.cpp

clean:
	rm *.o
