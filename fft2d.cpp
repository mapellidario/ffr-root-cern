/*
FFT in 2d
http://www.mathworks.it/it/help/images/fourier-transform.html
http://www.mathworks.it/it/help/matlab/ref/fftshift.html

g++ -o fft2d fft2d.cpp `root-config --cflags --glibs`
./fft2d 
*/

#include "TH2D.h"
#include "TVirtualFFT.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TROOT.h"

#include "TAxis.h"
#include "TH2.h"
#include "TArrayD.h"

//STL
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//funzioni per riscalare assi binnati
Double_t ScaleX2(Double_t x) ;
Double_t ScaleX1(Double_t x) ;
Double_t TraslateX(Double_t x);
Double_t ScaleY2(Double_t x) ;
Double_t ScaleY1(Double_t x) ;
Double_t TraslateY(Double_t y) ;
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t)) ;
void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t)) ;
void ScaleYaxis(TH1 *h, Double_t (*Scale)(Double_t)) ; //qui non usata

//funzione arbitraria da passare alla TF1, permettead esempio di avere funzion definite a tratti
double funzione_arbitraria (double *x, double *par) ;
int numero_parametri=0;

//parametri TH1 e TF1
Int_t xbins=400; //indica il numero di bin in cui voglio dividere l'asse x del mio istogramma
Int_t ybins=400;
Int_t n[2]={xbins,ybins}; //serve per la fft
double xmin=-2*TMath::Pi(); //valore minimo dell'assex degli isogrammi / funzioni
double xmax=2*TMath::Pi(); //valore massimo dell'assex degli isogrammi / funzioni
double ymin=-1*TMath::Pi(); //valore minimo dell'assey degli isogrammi / funzioni
double ymax=1*TMath::Pi(); //valore massimo dell'assey degli isogrammi / funzioni

//MAIN
int main (int numArg, char * listArg[]) {
	TApplication* theApp = new TApplication("App", &numArg, listArg);
	const bool BatchPrint = 1;
	if (BatchPrint) gROOT->SetBatch(); //Permette di editare immagini senza visualizzarle sullo schermo //STAMPA
	
//********* Histograms ********//
	//prepare the canvas for drawing
	TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 1000, 1000); //STAMPA
//	TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 2000, 2000);
//	myc->SetBatch();
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

	//-----------------------------------	
	//Represent the funcion you want to compute the fft of
	//A function to sample
//	TF2 *f1 = new TF2("f1", "sin(x)-sin(y)", xmin, xmax, ymin, ymax);
	TF2 *f1 = new TF2 ("f1", funzione_arbitraria, xmin, xmax, ymin, ymax, numero_parametri);
//	f1->Draw("contz"); //penso che si inventi un TGraph/TH1 su cui rappresentarla
	//representing the function through a histogram
	c1_1->cd();
	TH1 *h1 = new TH2D("h1", "f1 represented through a histogram", xbins, xmin, xmax, ybins, ymin, ymax);
	Double_t x,y;
	//Fill the histogram with function values
	for (Int_t i=1; i<=xbins; i++){
		for (Int_t j=1; j<=ybins; j++){
			x = (Double_t(i)/xbins)*(xmax-xmin)+xmin;
			y = (Double_t(j)/ybins)*(ymax-ymin)+ymin;
			h1->SetBinContent(i, j, f1->Eval(x,y));
		}		
	}
	/* if using a function with number of parameters different from 0.
	Double_t x[2], par[1];
	//Fill the histogram with function values
	for (Int_t i=1; i<=xbins; i++){
		for (Int_t j=1; j<=ybins; j++){
			x[0] = (Double_t(i)/xbins)*(xmax-xmin)+xmin;
			x[1] = (Double_t(j)/ybins)*(ymax-ymin)+ymin;
			par[0] = 1;
			h1->SetBinContent(i, j, f1->EvalPar(x,par));
		}		
	}
	*/
	h1->Draw("colz"); //disegna l'istogramma sullo stesso pad dove ho messo la funzione.
	//draw dirctly the TF2, but changing the resolution //Alternatively
	f1->SetNpx(1000);
	f1->SetNpy(1000);
//	f1->Draw(); 

	//-----------------------------------	
	//Compute the transform	
	TVirtualFFT::SetTransform(0);
	
	//look at the magnitude of the output
	c1_3->cd();
	c1_3->SetLogz();
	TH1 *hm =0;
	hm = h1->FFT(hm, "MAG P");
	hm->SetTitle("Magnitude of the 1st transform (with SetLogz())");
	hm->SetStats(kFALSE);
	hm->Draw("colz");
	//NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function (in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
	//di sicuro se min diverso da zero devo moltiplicare per (max-min)
	//non sono sicuro di ciò: se aumento il numero di bin, è normale che anche questo aumenti dello stesso fattore? INDAGARE
	//rescaling axis, as suggested above
	//rescaling a binned axis: using all the function listed before the main
	ScaleXaxis(hm, ScaleX1);
	ScaleYaxis(hm, ScaleY1);
	hm->Draw("colz");
	//rescaling the axis of the istogram which contains the information about the bins. aka multiplying the istogram
	hm->Scale(1./sqrt(xbins*ybins));
	hm->Draw("colz");
	
	//Look at the phase of the output	
	c1_4->cd();
	TH1 *hp = 0;
	hp = h1->FFT(hp, "PH P");
	hp->SetTitle("Phase of the 1st transform");
	hp->SetStats(kFALSE);
	hp->Draw("colz");
	//rescaling axis, but it was not suggested to. so i dont do it
	//rescaling a binned axis: using all the function listed before the main
	ScaleXaxis(hp, ScaleX1);
	ScaleYaxis(hp, ScaleY1);
	hp->Draw("colz");

	//-----------------------------------  
	//Use the following method to get the full output:
//	Double_t *re_full = new Double_t[xbins*ybins];
//	Double_t *im_full = new Double_t[xbins*ybins];
	Double_t *re_full = new Double_t[xbins*ybins];
	Double_t *im_full = new Double_t[xbins*ybins];	
	//That's the way to get the current transform object:
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	//Use the following method to get the full output:
	fft->GetPointsComplex(re_full,im_full);

	//-----------------------------------  
	//Now let's make a backward transform:
	c1_2->cd();
	TVirtualFFT *fft_back = TVirtualFFT::FFT(2, n, "C2R M K"); //non uso &n perchè ora n è già un puntatore
	fft_back->SetPointsComplex(re_full,im_full);
	fft_back->Transform();
	//Let's look at the output
	TH1 *hb = 0;
	hb = TH1::TransformHisto(fft_back,hb,"Re");
	hb->SetTitle("The backward transform result");
	hb->Draw("colz"); //non funziona benissimo: pecie la parte centrale è scialacquata
	//NOTE: here you get at the x-axes number of bins and not real values (in this case 25 bins has to be rescaled to a range between 0 and 4*Pi; also here the y-axes has to be rescaled (factor 1/bins) (VALE SOLO SE HO MIN=0)
	//per riscalare e avere min diverso da zero Bisogna usare: nuova_x = vecchia_x/n * (max-min)+min (implementato in ScaleX2)

	//rescaling axis, as suggested above
	//rescaling a binned axis: using all the function listed before the main
	ScaleXaxis(hb, ScaleX2);
	ScaleYaxis(hb, ScaleY2);
	//rescaling the axis of the istogram which contains the information about the bins. aka multiplyin the istogram
	hb->Scale(1./(xbins*ybins));
	hb->ResetStats(); //due to rescaling
	hb->SetStats(kFALSE);
	hb->Draw("colz");	

//----------------------------------------------------------------------------
	//reconstruct the magnitude hisogram to make it look like to the ones of matlab
	c1_5->cd();
	c1_5->SetLogz();

	TH1 *hfinal = new TH2D ("hfinal", "hfinal", hm->GetNbinsX(), hm->GetXaxis()->GetXmin(), hm->GetXaxis()->GetXmax(), 
																hm->GetNbinsY(), hm->GetYaxis()->GetXmin(), hm->GetYaxis()->GetXmax());
	for(int i=1; i<=xbins/2; i++){ //III quad -> I quad
		for(int j=1; j<=ybins/2; j++){
			hfinal->SetBinContent(i+xbins/2, j+ybins/2, hm->GetBinContent(i,j));
		}
	}	
	for(int i=xbins/2+1; i<=xbins; i++){ //I quad -> III quad
		for(int j=ybins/2+1; j<=ybins; j++){
			hfinal->SetBinContent(i-xbins/2, j-ybins/2, hm->GetBinContent(i,j));
		}
	}
	for(int i=xbins/2+1; i<=xbins; i++){ // IV quad -> II quad
		for(int j=1; j<=ybins/2; j++){
			hfinal->SetBinContent(i-xbins/2, j+ybins/2, hm->GetBinContent(i,j));
		}
	}
	for(int i=1; i<=xbins/2; i++){
		for(int j=ybins/2+1; j<=ybins; j++){ // II quad -> IV quad
			hfinal->SetBinContent(i+xbins/2, j-ybins/2, hm->GetBinContent(i,j));
		}
	}
	ScaleXaxis(hfinal, TraslateX); 
	ScaleYaxis(hfinal, TraslateY); 
	hfinal->SetTitle("Magnitude of the 1st transform (with SetLogz()) - Shifted (Matlab fftshift-like)");
	hfinal->SetStats(kFALSE);
	hfinal->Draw("colz");

	//reconstruct the phase hisogram to make it look like to the ones of matlab
	c1_6->cd();
	TH1 *hfinal2 = new TH2D ("hfinal2", "hfinal2", hp->GetNbinsX(), hp->GetXaxis()->GetXmin(), hp->GetXaxis()->GetXmax(), 
																hp->GetNbinsY(), hp->GetYaxis()->GetXmin(), hp->GetYaxis()->GetXmax());
	for(int i=1; i<=xbins/2; i++){ //III quad -> I quad
		for(int j=1; j<=ybins/2; j++){
			hfinal2->SetBinContent(i+xbins/2, j+ybins/2, hp->GetBinContent(i,j));
		}
	}	
	for(int i=xbins/2+1; i<=xbins; i++){ //I quad -> III quad
		for(int j=ybins/2+1; j<=ybins; j++){
			hfinal2->SetBinContent(i-xbins/2, j-ybins/2, hp->GetBinContent(i,j));
		}
	}
	for(int i=xbins/2+1; i<=xbins; i++){ // IV quad -> II quad
		for(int j=1; j<=ybins/2; j++){
			hfinal2->SetBinContent(i-xbins/2, j+ybins/2, hp->GetBinContent(i,j));
		}
	}
	for(int i=1; i<=xbins/2; i++){
		for(int j=ybins/2+1; j<=ybins; j++){ // II quad -> IV quad
			hfinal2->SetBinContent(i+xbins/2, j-ybins/2, hp->GetBinContent(i,j));
		}
	}
	ScaleXaxis(hfinal2, TraslateX); 
	ScaleYaxis(hfinal2, TraslateY); 
	hfinal2->SetTitle("Phase of the 1st transform - Shifted (Matlab fftshift-like)");
	hfinal2->SetStats(kFALSE);
	hfinal2->Draw("colz");

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
		in[i] =  f1->Eval(x);
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
	
	delete [] in;
*/

	delete [] re_full;
	delete [] im_full;

	//Stampa Canvas //STAMPA
	if (BatchPrint) {
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
/*		std::string title_tex;  
		std::string C = ".tex";
		title_tex = aux + C ; //creo il nome da assegnare al file .tex
		myc->Print(title_tex.c_str());
*/
	}
	if(!BatchPrint) theApp->Run(); //STAMPA
	return 0;
}

//Funzioni
	Double_t ScaleX2(Double_t x)	{
		Double_t v;
		v = x/xbins * (xmax-xmin)+xmin ; 
		return v;
	}
	
	Double_t ScaleX1(Double_t x)	{
		Double_t v;
		v = x / (xmax-xmin) ; 
		return v;
	}
	
	Double_t TraslateX(Double_t x)	{
		Double_t v;
		v = x - (xbins/2) / (xmax-xmin) ; 
		return v;
	}
	
	Double_t ScaleY2(Double_t y)	{
		Double_t v;
		v = y/ybins * (ymax-ymin)+ymin ; 
		return v;
	}
	
	Double_t ScaleY1(Double_t y)	{
		Double_t v;
		v = y / (ymax-ymin) ; 
		return v;
	}	
	
	Double_t TraslateY(Double_t y)	{
		Double_t v;
		v = y - (ybins/2) / (ymax-ymin) ; 
		return v;
	}

	void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))	{
		if (!a) return; // just a precaution
		if (a->GetXbins()->GetSize())		{
			// an axis with variable bins
			// note: bins must remain in increasing order, hence the "Scale"
			// function must be strictly (monotonically) increasing
			TArrayD X(*(a->GetXbins()));
			for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
			a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
		}
		else		{
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

	void ScaleYaxis(TH1 *h, Double_t (*Scale)(Double_t))	{
		if (!h) return; // just a precaution
		ScaleAxis(h->GetYaxis(), Scale);
		return;
	}
 
	//--------------------------------------------
double funzione_arbitraria (double *x, double *par) {
	if (x[0]>-0.3 && x[0]<0.3 && x[1]>x[0]/6.-1./10. && x[1]<x[0]/6.+1./10.) return 1. ;
	else return 0.;
}
