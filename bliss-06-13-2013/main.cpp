/*

Copyright or © or Copr. Leo Ducas and Tancrede Lepoint.

Leo.Ducas@ens.fr and Tancrede.Lepoint@ens.fr

This software is a computer program whose purpose is to provide to the 
research community a proof-of-concept implementation of the BLISS 
digital signature scheme of Ducas, Durmus, Lepoint and Lyubashevsky 
appeared at Crypto 2013.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/

#include <iostream>
#include "Setup.h"
#include "KeyGen.h"
#include "Sampler.h"
#include "Sign.h"
#include "Verify.h"
#include "Entropy.h"
#include <ctime>
#include <mpfr.h>
#include <fstream>
#include "NTL/ZZ_pE.h"


// Set to 1 if you want to verify
#define VERIFY 1

// Set to 1 if you want to print an array of signatures in SAGE format
#define PRINT 0


int main() {

	long i;
    time_t timeStart;
    timespec start, finish;
    double cpu_time, monotonic_time;

#define do_time_start()	do { \
    timeStart = clock(); \
    clock_gettime(CLOCK_MONOTONIC, &start); } while(0)
#define do_time_finish(m) do { \
    cpu_time = (static_cast<double>(clock()-timeStart))/CLOCKS_PER_SEC/(m); \
    clock_gettime(CLOCK_MONOTONIC, &finish); \
    monotonic_time = (finish.tv_sec - start.tv_sec)/(m) + \
        (finish.tv_nsec - start.tv_nsec) / 1000000000.0 / (m); } while(0)
#define elapsed	cpu_time << " s (CPU) / " << monotonic_time << " s (wall)"

    /* Comparison randomness 
    Entropy random;

    do_time_start();
    for (i=0; i<1000000; i++)
    	random.getURandomLong();
    do_time_finish(1.0);
	std::cout << "Urandom: " << elapsed << std::endl;

	do_time_start();
    for (i=0; i<1000000; i++)
    	random.getRandomLong();
    do_time_finish(1.0);
	std::cout << "Sha512 : " << elapsed << std::endl;
	//*/

	/* Check randomness
	Entropy random;

	long histochar[256];
	for (i=0; i<256; i++) histochar[i]=0;
	for (i=0; i<5000000; i++)
		histochar[random.getRandomChar()]++;

	for (i=0; i<256; i++) std::cout << i << ": "<< histochar[i] << std::endl;
	
	//*/

	/* Setup */
    Setup setup;
	Entropy random;
	Sampler sampler(sigma, alpha_rejection, &random);

	//*/


	/* Test Gaussian Sampling

	Entropy random;
	Sampler sampler(81, .5, &random);
	do_time_start()	;

	for (int i=0; i<5000000; i++){
	  sampler.SamplerGaussian();
	}
	do_time_finish(1.0);
	std::cout << "Sampling 100 000 entries : " << elapsed << std::endl;

	//*/



	std::cout << "CLASS-" << CLASS << std::endl;

	/* KeyGen */ 
	do_time_start();
	KeyGen key(setup, &random);
	/*KeyGen* kkey;
	for (i=0; i<20; i++) kkey = new KeyGen(setup, &random);
	*/
	do_time_finish(1.0);
	std::cout << "KeyGen: " << elapsed << std::endl;
	//*/

	/* Verification AS = T */ 

	NTL::ZZ_p::init(NTL::to_ZZ(2*Q));
	NTL::ZZ_pX phi;
	conv(phi, setup.get_phi());
	NTL::ZZ_pE::init(phi);

	NTL::ZZ_pE a1, a2, s1, s2;
	NTL::ZZ_pX pa1, pa2, ps1, ps2;
	for (i=0; i<N; i++) {
		NTL::SetCoeff(pa1, i, key.pk.a1[key.pk.offset+i]);
	}
	NTL::SetCoeff(pa2, 0, key.pk.a2);
	conv(ps1, key.sk.s1);
	conv(ps2, key.sk.s2);
	conv(a1, pa1);
	conv(a2, pa2);
	conv(s1, ps1);
	conv(s2, ps2);

	std::cout << "a1=" << a1 << std::endl;
	std::cout << "a2=" << a2 << std::endl;
	std::cout << "s1=" << s1 << std::endl;
	std::cout << "s2=" << s2 << std::endl;

	std::cout << "a*s=" << (a1*s1+a2*s2) << std::endl;
	//*/

	/* Sign */
	Sign sign(setup, &sampler, &random);
	Verify verify;

	std::string message = "Lorem ipsum dolor sit amet, consectetur adipiscing \
	elit. Nulla nec posuere mauris, eget facilisis nunc. Maecenas id libero at \
	tortor luctus volutpat. Nulla consequat tellus non ligula tincidunt, sed \
	tristique lectus volutpat. Curabitur vestibulum pretium blandit. Nunc nec \
	risus molestie urna elementum pellentesque. Nam id facilisis est, a \
	ultricies metus. Vestibulum sit amet nibh nec justo varius convallis id \
	nec justo. Maecenas id risus dignissim, volutpat nisi nec, fringilla elit. \
	Donec varius sem id risus sagittis, eu aliquam sem feugiat. Nulla non leo \
	tempus, pulvinar velit ut, euismod dui. Aliquam leo nunc, auctor et nulla \
	quis, dapibus vulputate quam. Morbi ut faucibus nisl, vitae vehicula elit. \
	Ut sed purus enim. Maecenas sit amet mi volutpat, rutrum velit quis, \
	volutpat quam. Vestibulum lacus tellus, placerat vitae magna vitae, \
	laoreet feugiat magna. Aliquam nulla lacus, scelerisque eu mauris a, \
	viverra suscipit ligula.";
	
	long totalTests = 0;
	long nbSign = 10000;
	//*/

	#if VERIFY
	timespec verifystart, verifyfinish;
	double verifytime = 0.0;
	int nb_verif = 0;
	#endif

#if PRINT
	std::cerr << "CLASS-"<< CLASS << " Dim : "  << N << std::endl;

std::cout << "L=[";
#endif
do_time_start();
for (i=0; i<nbSign; i++) {
totalTests += sign.signMessage(key.pk, key.sk, message);
#if PRINT
long k;
std::cout << "(";

std::cout << "[";
for (k=0; k<N; k++)
{
std::cout << sign.signOutput.z1[k];
if (k<N-1)
std::cout << ",";
}
std::cout << "]";

std::cout << ",[";
for (k=0; k<N; k++)
{
std::cout << sign.signOutput.z2Carry[k];


if (k<N-1)
std::cout << ",";
}
std::cout << "]";

std::cout << ",[";
for (k=0; k<N; k++)
{
long l;
bool out=false;
for (l=0; l<kappa; l++)
{
if (sign.signOutput.indicesC[l]==k)
{
out=true;
std::cout << "1";
}
}
if (!out)
std::cout << "0";
if (k<N-1)
std::cout << ",";
}
std::cout << "])";

if (i<nbSign-1)
std::cout << ",";
#endif
#if VERIFY
clock_gettime(CLOCK_MONOTONIC, &verifystart);
nb_verif += verify.verifyMessage(key.pk, sign.signOutput, message);

clock_gettime(CLOCK_MONOTONIC, &verifyfinish);
verifytime += (verifyfinish.tv_sec - verifystart.tv_sec) + (verifyfinish.tv_nsec - verifystart.tv_nsec) / 1000000000.0;
#endif
}
#if PRINT
std::cout << "]" << std::endl;
#endif

	do_time_finish((double)nbSign);
	std::cout << "Sign: " << elapsed << std::endl;
	#if VERIFY
	std::cout << "Including Verify: " << verifytime/nbSign << " s (wall) ; Number of verified signatures: " << nb_verif << "/" << nbSign << std::endl;
	#endif

	//*/
	
	//for (i=0;i<100;i++)
	//  std::cerr << random.getRandomLong() << "\t";

	/* Sampler Gaussian 
	long histo2[20000];
	for (i=0; i<20000; i++)
		histo2[i] = 0;
	do_time_start();
	for (i=0; i<50000000; i++) {
		histo2[10000+sampler.SamplerGaussian()]++;
	}
	do_time_finish(1000.0);
	std::cout << "SamplerGaussian: " << elapsed << std::endl;
	
	std::ofstream myfile;
	myfile.open ("histo.dat");

	for (i=-999; i<1000; i++)
		myfile << i << "\t" << histo2[10000+i] << std::endl;
	myfile.close();
	
	double mean=0,
	  stddev=0,
	  NN=0;
	for (int i=-10000; i<10000; i++)
	  {
	    mean += (double)i*histo2[10000+i];
	    NN+=histo2[10000+i];
	  }
	mean /= NN;
	
	std::cout << "Mean = " << mean << std::endl;
	
	for (i=-10000; i<10000; i++)
	{
		stddev += (double) (i-mean)*(i-mean)*histo2[10000+i];	
	}
	stddev /= NN;
	std::cout << "StdDev^2=" << stddev << std::endl;
	//*/

	/* Sampler Positive Binary 
	long histo[11];
	for (i=0; i<11; i++)
		histo[i] = 0;
	
	do_time_start();
	for (i=0; i<10000000; i++)
		histo[sampler.SamplerPosBin()]++;
	do_time_finish(10000.0);
	std::cout << "SamplerPosBin: " << elapsed << std::endl;
	for (i=0; i<11; i++)
		std::cout << histo[i]<<std::endl;
	//*/

	/* SamplerBernoulli(x) 
	unsigned char x[16];
	x[0]=64;
	for (i=1; i<16; i++) x[i] = 0;
	
	long histo3[2];
	histo3[0]=0;
	histo3[1]=0;
	
	for (i=0; i<100000; i++)
		histo3[sampler.SamplerBer(x)]++;
	
	std::cout << histo3[0] << std::endl;
	std::cout << histo3[1] << std::endl;
	//*/

	/*	Sampler Bernoulli Exponential 
	long histo4[2];
	histo4[0]=0; histo4[1]=0;
	
	for (i=0; i<1000000; i++)
		histo4[sampler.SamplerBerExp(7500)]++;
	
	std::cout << histo4[0] << std::endl;
	std::cout << histo4[1] << std::endl;
	//*/

	/*	Sampler Bernoulli ExpM 
	long histo5[2];
	histo5[0]=0; histo5[1]=0;
	
	for (i=0; i<1000000; i++)
		histo5[sampler.SamplerBerExpM(1000)]++;//histo5[sampler.SamplerBerExpM(10000)]++;
	
	std::cout << histo5[0] << std::endl;
	std::cout << histo5[1] << std::endl;
	//*/

	/*	Sampler Bernoulli Cosh 
	long histo6[2];
	histo6[0]=0; histo6[1]=0;
	
	for (i=0; i<1000000; i++)
		histo6[sampler.SamplerBerCosh(10000)]++;
	
	std::cout << histo6[0] << std::endl;
	std::cout << histo6[1] << std::endl;
	//*/
    return 0;
}

