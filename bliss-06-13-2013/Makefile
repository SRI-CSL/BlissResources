# Copyright or © or Copr. Leo Ducas and Tancrede Lepoint.
#
# Leo.Ducas@ens.fr and Tancrede.Lepoint@ens.fr
# 
# This software is a computer program whose purpose is to provide to the 
# research community a proof-of-concept implementation of the BLISS 
# digital signature scheme of Ducas, Durmus, Lepoint and Lyubashevsky 
# appeared at Crypto 2013.
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

CXX	= g++
IFLAGS = 
LDFLAGS = -lgmpxx -lgmp -lrt -lntl -lmpfr -lgmp -lcrypto
CXXFLAGS= -march=native -mtune=native -O3

all: Makefile Entropy.cpp Entropy.h Fft.cpp Fft.h KeyGen.cpp KeyGen.h main.cpp Parameters.h Sampler.cpp Sampler.h Setup.cpp Setup.h Sign.cpp Sign.h Verify.cpp Verify.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o Entropy.o Entropy.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o Fft.o Fft.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o KeyGen.o KeyGen.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o main.o main.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o Sampler.o Sampler.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o Setup.o Setup.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o Sign.o Sign.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o Verify.o Verify.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o bliss main.o Entropy.o Fft.o KeyGen.o Sampler.o Setup.o Sign.o Verify.o $(LDFLAGS) 

clean:
	rm -f *.o bliss histo.dat 
