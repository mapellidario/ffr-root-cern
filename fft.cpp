/*
//This tutorial illustrates the Fast Fourier Transforms interface in ROOT.
//FFT transform types provided in ROOT:
// - "C2CFORWARD" - a complex input/output discrete Fourier transform (DFT) 
//						in one or more dimensions, -1 in the exponent
// - "C2CBACKWARD"- a complex input/output discrete Fourier transform (DFT) 
//						in one or more dimensions, +1 in the exponent
// - "R2C"		- a real-input/complex-output discrete Fourier transform (DFT)
//						in one or more dimensions,
// - "C2R"		- inverse transforms to "R2C", taking complex input 
//						(storing the non-redundant half of a logically Hermitian array) 
//						to real output
// - "R2HC"		- a real-input DFT with output in Â¡ÃˆhalfcomplexÂ¡Ã‰ format, 
//						i.e. real and imaginary parts for a transform of size n stored as
//						r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
// - "HC2R"		- computes the reverse of FFTW_R2HC, above
// - "DHT"		- computes a discrete Hartley transform
// Sine/cosine transforms:
//  DCT-I  (REDFT00 in FFTW3 notation)
//  DCT-II (REDFT10 in FFTW3 notation)
//  DCT-III(REDFT01 in FFTW3 notation)
//  DCT-IV (REDFT11 in FFTW3 notation)
//  DST-I  (RODFT00 in FFTW3 notation)
//  DST-II (RODFT10 in FFTW3 notation)
//  DST-III(RODFT01 in FFTW3 notation)
//  DST-IV (RODFT11 in FFTW3 notation)
//First part of the tutorial shows how to transform the histograms
//Second part shows how to transform the data arrays directly
//Authors: Anna Kreshuk and Jens Hoffmann
Elaborato un pochino: riscalato gli assi (binnati e non) eriordinato l'ordine dei comandi

fftw oficial manual: http://www.fftw.org/fftw3_doc/Multi_002ddimensional-Transforms.html#Multi_002ddimensional-Transforms
The multi-dimensional transforms of FFTW, in general, compute simply the separable product of the given 1d transform along each dimension of the array. Since each of these transforms is unnormalized, computing the forward followed by the backward/inverse multi-dimensional transform will result in the original array scaled by the product of the normalization factors for each dimension (e.g. the product of the dimension sizes, for a multi-dimensional DFT). 

g++ -o fft fft.cpp `root-config --cflags --glibs`
./fft 
*/

#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TROOT.h"

#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"

//STL
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//funzioni per riscalare assi binnati
Double_t ScaleX2(Double_t x) ;
Double_t ScaleX1(Double_t x) ;
Double_t TraslateX(Double_t x) ;
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t)) ;
void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t)) ;

//funzione arbitraria da passare alla TF1, permettead esempio di avere funzion definite a tratti
double funzione_arbitraria (double *x, double *par) ;
int numero_parametri=0;

//parametri TH1 e TF1
Int_t n=250; //indica il numero di bin (in realtà i bin sono n+1) in cui voglio dividere l'asse x del mio istogramma
double min=-4*TMath::Pi(); //valore minimo dell'assex degli isogrammi / funzioni
double max=4*TMath::Pi(); //valore massimo dell'assex degli isogrammi / funzioni


int main (int numArg, char * listArg[])
{
	TApplication* theApp = new TApplication("App", &numArg, listArg);
//	gROOT->SetBatch(); //Permette di editare immagini senza visualizzarle sullo schermo
	
//********* Histograms ********//
	//prepare the canvas for drawing
	TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 1000, 1000);
//	TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 2000, 2000); //STAMPA
	myc->SetFillColor(45);
	TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.67,0.49,0.99);
	TPad *c1_2 = new TPad("c1_2", "c1_2",0.51,0.67,0.99,0.99);
	TPad *c1_3 = new TPad("c1_3", "c1_3",0.01,0.34,0.49,0.65);
	TPad *c1_4 = new TPad("c1_4", "c1_4",0.51,0.34,0.99,0.65);
	TPad *c1_5 = new TPad("c1_5", "c1_5",0.01,0.01,0.49,0.32);
	TPad *c1_6 = new TPad("c1_6", "c1_6",0.51,0.01,0.99,0.32);
	c1_1->Draw();
	c1_2->Draw();
	c1_3->Draw();
	c1_4->Draw();
	c1_5->Draw();
	c1_6->Draw();
	c1_1->SetFillColor(30);
	c1_1->SetFrameFillColor(42);
	c1_2->SetFillColor(30);
	c1_2->SetFrameFillColor(42);
	c1_3->SetFillColor(30);
	c1_3->SetFrameFillColor(42);
	c1_4->SetFillColor(30);
	c1_4->SetFrameFillColor(42);
	c1_5->SetFillColor(30);
	c1_5->SetFrameFillColor(42);
	c1_6->SetFillColor(30);
	c1_6->SetFrameFillColor(42);

	//-----------------------------------	
	//Represent the funcion you want to compute the fft of
	c1_1->cd();
	TH1::AddDirectory(kFALSE);	
	//A function to sample
	TF1 *fsin = new TF1("fsin", "sin(x)+sin(2*x)+sin(0.5*x)+1", min, max);
//	TF1 *fsin = new TF1("fsin", "sin(0.5*x)+3*sin(x)+sin(2*x)", min, max);
//	TF1 *fsin = new TF1 ("fsin", funzione_arbitraria, min, max, numero_parametri);
	fsin->Draw(); //penso che si inventi un TGraph/TH1 su cui rappresentarla

	TH1D *hsin = new TH1D("hsin", "hsin", n, min, max);//NB rispetto a quello originale di ROOT ho cambiato il numero di bin
	Double_t x;
	//Fill the histogram with function values
	for (Int_t i=1; i<=n; i++){ //NB rispetto a quello originale di ROOT ho cambiato gli estremi del parametro i
		x = (Double_t(i)/n)*(max-min)+min; //ora posso avere min diverso da zero
		hsin->SetBinContent(i, fsin->Eval(x));
	}
	hsin->Draw("same"); //disegna l'istogramma sullo stesso pad dove ho messo la funzione.
	fsin->GetXaxis()->SetLabelSize(0.05);
	fsin->GetYaxis()->SetLabelSize(0.05);

	//-----------------------------------	
	//Compute the transform	
	TVirtualFFT::SetTransform(0);
	
	//look at the magnitude of the output
	c1_3->cd();
	c1_3->SetLogz();
	TH1 *hm =0;
	hm = hsin->FFT(hm, "MAG");
	hm->SetTitle("Magnitude of the 1st transform");
	hm->SetStats(kFALSE);
	hm->Draw();
	//NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function (in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
	//di sicuro se min diverso da zero devo moltiplicare per (max-min)
	//rescaling axis, as suggested above
	//rescaling a binned axis: using all the function listed before the main
	ScaleXaxis(hm, ScaleX1);
	hm->Draw();
	//rescaling the axis of the istogram which contains the information about the bins. aka multiplying the istogram
	hm->Scale(1./sqrt(n));
	hm->Draw();	
	hm->ResetStats(); //due to rescaling
	
	//Look at the phase of the output	
	c1_4->cd();
	TH1 *hp = 0;
	hp = hsin->FFT(hp, "PH");
	hp->SetStats(kFALSE);
	hp->SetTitle("Phase of the 1st transform");
	hp->Draw();
	//rescaling axis, even if it was not suggested to.
	//rescaling a binned axis: using all the function listed before the main
	ScaleXaxis(hp, ScaleX1);
	hp->Draw();
	hp->ResetStats(); //due to rescaling
 
	//-----------------------------------  
	//OUTPUT
	//That's the way to get the current transform object:
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	//Use the following method to get the full output:
	Double_t *re_full = new Double_t[n];
	Double_t *im_full = new Double_t[n];
	fft->GetPointsComplex(re_full,im_full);
	
	//Look at the DC component and the Nyquist harmonic:  
	//Use the following method to get just one point of the output
	Double_t re, im;
	fft->GetPointComplex(0, re, im);
	printf("1st transform: DC component: %f\n", re);
	fft->GetPointComplex(n/2+1, re, im);
	printf("1st transform: Nyquist harmonic: %f\n", re);



	//-----------------------------------  
	//Now let's make a backward transform:
	c1_2->cd();
	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
	fft_back->SetPointsComplex(re_full,im_full);
	fft_back->Transform();

	//Let's look at the output
	TH1 *hb = 0;
	hb = TH1::TransformHisto(fft_back,hb,"Re");
	hb->SetTitle("The backward transform result");
  	hb->SetStats(kFALSE);
	hb->Draw();
	//NOTE: here you get at the x-axes number of bins and not real values (in this case 25 bins has to be rescaled to a range between 0 and 4*Pi; also here the y-axes has to be rescaled (factor 1/bins) (VALE SOLO SE HO MIN=0)
	//per riscalare e avere min diverso da zero Bisogna usare: nuova_x = vecchia_x/n * (max-min)+min (implementato in ScaleX2)

	//rescaling axis, as suggested above
	//rescaling a binned axis: using all the function listed before the main
	ScaleXaxis(hb, ScaleX2);
	hb->Draw();
	//rescaling the axis of the istogram which contains the information about the bins. aka multiplyin the istogram
	hb->Scale(1./n);
	hb->Draw();	
	hb->ResetStats(); //due to rescaling

	//superimpose the original function to the backward transform, checking if it worked well
	fsin->Draw("same"); //penso che si inventi un TGraph/TH1 su cui rappresentarla

	//--------------------------------------------------------------------------------
	//reconstruct the magnitude hisogram to make it look like to the ones of matlab
	//huston: we got a problem
	c1_5->cd();
	c1_5->SetLogz();

	TH1 *hfinal = new TH1D ("hfinal", "hfinal", hm->GetNbinsX(), hm->GetXaxis()->GetXmin(), hm->GetXaxis()->GetXmax());
	for(int i=1; i<=n/2; i++){ //III quad -> I quad
		hfinal->SetBinContent(i+n/2, hm->GetBinContent(i));
	}	
	for(int i=n/2+1; i<=n; i++){ //I quad -> III quad
		hfinal->SetBinContent(i-n/2, hm->GetBinContent(i));
	}
	ScaleXaxis(hfinal, TraslateX); //INDINSPENSABILE
	hfinal->SetTitle("Magnitude of the 1st transform (with SetLogz()) - Shifted (Matlab fftshift-like)");
	hfinal->SetStats(kFALSE);
	hfinal->Draw("colz");

	//reconstruct the phase hisogram to make it look like to the ones of matlab
	//huston: we got a problem
	c1_6->cd();

	TH1 *hfinal1 = new TH1D ("hfinal1", "hfinal1", hp->GetNbinsX(), hp->GetXaxis()->GetXmin(), hp->GetXaxis()->GetXmax());
	for(int i=1; i<=n/2; i++){ 
		hfinal1->SetBinContent(i+n/2, hp->GetBinContent(i));
	}	
	for(int i=n/2+1; i<=n; i++){
		hfinal1->SetBinContent(i-n/2, hp->GetBinContent(i));
	}
	ScaleXaxis(hfinal1, TraslateX); //INDINSPENSABILE
	hfinal1->SetTitle("Phase of the 1st transform - Shifted (Matlab fftshift-like)");
	hfinal1->SetStats(kFALSE);
	hfinal1->Draw("colz");

	delete fft_back;
	fft_back=0;

//********* Data array - same transform ********//
/*  //Allocate an array big enough to hold the transform output
	//Transform output in 1d contains, for a transform of size N, 
	//N/2+1 complex numbers, i.e. 2*(N/2+1) real numbers
	//our transform is of size n+1, because the histogram has n+1 bins

	Double_t *in = new Double_t[2*((n+1)/2+1)];
	Double_t re_2,im_2;
	for (Int_t i=0; i<=n; i++){
		x = (Double_t(i)/n)*(max-min)+min; //ora posso avere min diverso da zero
		in[i] =  fsin->Eval(x);
	}

	//Make our own TVirtualFFT object (using option "K")
	//Third parameter (option) consists of 3 parts:
	//-transform type:
	// real input/complex output in our case
	//-transform flag: 
	// the amount of time spent in planning
	// the transform (see TVirtualFFT class description)
	//-to create a new TVirtualFFT object (option "K") or use the global (default)
	Int_t n_size = n+1;
	TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
	if (!fft_own) return 0;
	fft_own->SetPoints(in);
	fft_own->Transform();

	//Copy all the output points:
	fft_own->GetPoints(in);
	//Draw the real part of the output
	c1_5->cd();
	TH1 *hr = 0;
	hr = TH1::TransformHisto(fft_own, hr, "RE");
	hr->SetTitle("Real part of the 3rd (array) tranfsorm");
	hr->Draw();
	hr->SetStats(kFALSE);
	hr->GetXaxis()->SetLabelSize(0.05);
	hr->GetYaxis()->SetLabelSize(0.05);
	c1_6->cd();
	TH1 *him = 0;
	him = TH1::TransformHisto(fft_own, him, "IM");
	him->SetTitle("Im. part of the 3rd (array) transform");
	him->Draw();
	him->SetStats(kFALSE);
	him->GetXaxis()->SetLabelSize(0.05);
	him->GetYaxis()->SetLabelSize(0.05);

	//Now let's make another transform of the same size
	myc->cd();
	//The same transform object can be used, as the size and the type of the transform haven't changed
	TF1 *fcos = new TF1("fcos", "cos(x)+cos(0.5*x)+cos(2*x)+1", min, max);
	for (Int_t i=0; i<=n; i++){
		x = (Double_t(i)/n)*(4*TMath::Pi());
		in[i] =  fcos->Eval(x);
	}
	fft_own->SetPoints(in);
	fft_own->Transform();
	fft_own->GetPointComplex(0, re_2, im_2);
	printf("2nd transform: DC component: %f\n", re_2);
	fft_own->GetPointComplex(n/2+1, re_2, im_2);
	printf("2nd transform: Nyquist harmonic: %f\n", re_2);
	delete fft_own;
	delete [] in;
	delete [] re_full;
	delete [] im_full;
*/	
	
	
	//Stampa Canvas
	std::string aux = listArg[0] ;
	myc->Modified();
	myc->Update();
	myc->Draw();
	//decommentare per generare il file .png
	std::string title_png;
	std::string B = ".png";
	title_png = aux + B; //creo il nome da assegnare al file .png
	myc->Print(title_png.c_str(),"png");
	std::string title_pdf;
	std::string C = ".pdf";
	title_pdf = aux + C; //creo il nome da assegnare al file .png
	myc->Print(title_pdf.c_str(),"pdf");
	//non genero il file .tex, poichè pesa 80 MB
/*	std::string title_tex;  
	std::string C = ".tex";
	title_tex = aux + C ; //creo il nome da assegnare al file .tex
	myc->Print(title_tex.c_str());
*/	
	theApp->Run(); //TApplication
	return 0;
}

//Funzioni
	Double_t ScaleX2(Double_t x)	{
		Double_t v;
		v = x/n * (max-min)+min ; // posso riscalare e avere min diverso da zero
		return v;
	}
	
	Double_t ScaleX1(Double_t x)	{
		Double_t v;
		v = x / (max-min) ; // "linear scaling" function example
		return v;
	}
	
	Double_t TraslateX(Double_t x)	{
		Double_t v;
		v = x - (n/2) / (max-min) ; 
		return v;
	}

	void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))	{
		if (!a) return; // just a precaution
		if (a->GetXbins()->GetSize())	{
			// an axis with variable bins
			// note: bins must remain in increasing order, hence the "Scale"
			// function must be strictly (monotonically) increasing
			TArrayD X(*(a->GetXbins()));
			for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
			a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
		}
		else	{
			// an axis with fix bins
			// note: we modify Xmin and Xmax only, hence the "Scale" function
			// must be linear (and Xmax must remain greater than Xmin)
			a->Set( a->GetNbins(),
						Scale(a->GetXmin()), // new Xmin
						Scale(a->GetXmax()) ); // new Xmax
		}
		return;
	}

	void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t))	{
		if (!h) return; // just a precaution
		ScaleAxis(h->GetXaxis(), Scale);
		return;
	}
	
	//--------------------------------------------
double funzione_arbitraria (double *x, double *par) {
	if (x[0]<-1. || x[0]>2.) return 0. ;
	else return 1.;
}

