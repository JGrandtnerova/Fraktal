#include <iostream>
#include <fstream> /* pre ofstream */
#include <complex> /* for fractal arithmetic */
#include <omp.h>

typedef std::complex<double> COMPLEX;

class linear2d_function {
public:
	double a,b,c;
	void set(double a_,double b_,double c_) {a=a_;b=b_;c=c_;}
	linear2d_function(double a_,double b_,double c_) {set(a_,b_,c_);}
	double evaluate(double x,double y) const {return x*a+y*b+c;}
};

class pixel { /* onscreen RGB pixel */
public:
	unsigned char r,g,b;
	pixel() {}
	pixel(const COMPLEX &z) {
		r=(unsigned char)(z.real()*(256/2.0));
		g=(unsigned char)(z.imag()*(256/2.0));
		b=(unsigned char)(((z.real()*z.real()+z.imag()*z.imag()))*256);
	}
};


int main(void)
{
	// velkost obrazka
	int wid=1350, ht=1256;

	// Create a PPM output image file header:
	std::ofstream out("vystup.ppm",std::ios_base::binary);
	out<<"P6\n"
	   <<wid<<" "<<ht<<"\n"
	   <<"255\n";

    int tid; //premenna na pocitanie aktual. threadu - z cviceni
    double start_time = omp_get_wtime(); //nastavenie start time pocitadla casu
    int nthreads; //premenna pre zvolenie poctu threadov uzivatelom
    int max_count; // premenna na zvolenie poctu iteracii uzivatelom

    printf("zadaj pocet iteracii \n");
    scanf("%d",&max_count);
    printf("zadaj pocet threadov \n");
    scanf("%d",&nthreads);

    omp_set_num_threads(nthreads);

	// Set up coordinate system to render the Mandelbrot Set:
	double scale=3.0/wid;
	linear2d_function fx(scale,0.0,-scale*wid/2); // returns c given pixels
	linear2d_function fy(0.0,scale,-1.0);
	pixel *output=new pixel[wid*ht];


  /* Fork a team of threads - z cviceni*/
#pragma omp parallel private(tid) // nastavenie teamu threadov (a mastra), premennu tid vyuzivam pri zobrazeni aktualne beziaceho vlakna
  {
   tid = omp_get_thread_num();
    /* Only master thread does this - z cviceni, zostalo to na debug :)*/
    if (tid == 0)
      {
	nthreads = omp_get_num_threads();
	printf("Number of threads = %d\n", nthreads);
      }

#pragma omp for // tu bezia jednotlive thready pri pocitani vykreslovania
	for (int y=0;y<ht;y++) {
    //printf("Som thread c: %d \n",tid); debug
        for (int x=0;x<wid;x++) {
        //printf("Som thread c: %d \n",tid); debug

        //nastavenie pixelov pre vykreslenie
            COMPLEX c(fx.evaluate(x,y),fy.evaluate(x,y));
            COMPLEX z(0.0);
            int count;
            for (count=0;count<max_count;count++) { // vyuzitie premennej zvolenej uzivatelom na pocet iteracii v tomto for cykle
                z=z*z+c;
                if ((z.real()*z.real()+z.imag()*z.imag())>=4.0) break;
            }

        // Nastavenie farby pixelu do outputu
            unsigned char r,g,b;
            r=(unsigned char)(z.real()*(256/2.0));
            g=(unsigned char)(z.imag()*(256/2.0));
            b=(unsigned char)(((z.real()*z.real()+z.imag()*z.imag()))*256);
            output[x+y*wid]=pixel(z);
        }
    };
    printf("Som thread c: %d \n",tid); // kvoli prehladnosti ci fakt ide tolko threadov ako ma
  }
  /* End of parallel region */

	out.write((char *)&output[0],sizeof(pixel)*wid*ht); // output - vykreslenie obrazka


    double time = omp_get_wtime() - start_time; // vypocet casu behu programu- startovaci minus aktualny

    printf("Time: %lf\n",time); //vypis casu behu programu
	return 0;
}
