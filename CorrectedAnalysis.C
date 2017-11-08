#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <math.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TStyle.h"
#include "TColor.h"
#include <sstream> 
#include <fstream>
#include <stdlib.h>
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TLatex.h"
using namespace std;

	Double_t conc;
	TString name;
	vector<Int_t> cathode_number;
	vector<TString> cathode_name;
	vector<Double_t> cathode_concentration;
	
void CorrectedAnalysis()
{
	vector<vector<Double_t>> standard_list(100);
	vector<TString> standard_name_list;
	vector<vector<Double_t>> run_list(100);
	vector<TString> run_name_list;
	vector<Double_t> concentration_list;
	vector<Int_t> run_files;
	vector<Int_t> cathode_files;

	Int_t numtemp;
	TString nametemp;
	Double_t conctemp;
	Int_t runnumtemp;
	Int_t cathodetemp;
	Int_t cathodetemp2=-1;
	Int_t n=-1;
	Int_t j=-1;
	Int_t totalruns;
	Int_t totalstans;
	
	std::ifstream cathodefile;
	std::ifstream runfile;
	cathodefile.open("CathodeList.txt");
	runfile.open("RunList.txt");
	while(1)
	{
		cathodefile >> numtemp >> nametemp >> conctemp;
		cathode_number.push_back(numtemp);
		cathode_name.push_back(nametemp);
		cathode_concentration.push_back(conctemp);
		if(!cathodefile.good()) break;
	}
	while(1)
	{
		runfile >> runnumtemp >> cathodetemp;
		run_files.push_back(runnumtemp);
		cathode_files.push_back(cathodetemp);
		if(!runfile.good()) break;
	}
	cathodefile.close();
	runfile.close();
	for(int q=0;q<run_files.size();q++)
	{
		cathodetemp=cathode_files[q];
		runnumtemp=run_files[q];
		GetCathode(cathodetemp);
		if(name == "Standard")
		{
		if(cathodetemp != cathodetemp2)
		{
			n++;
			standard_name_list.push_back(name);
			concentration_list.push_back(conc);
		}
			totalstans++;
			standard_list[n].push_back(runnumtemp);
			cathodetemp2 = cathodetemp;
		}
		else if(name != "Standard")
		{
		if(cathodetemp != cathodetemp2)
		{
			j++;
			run_name_list.push_back(name);
		}
			totalruns++;
			run_list[j].push_back(runnumtemp);
			cathodetemp2 = cathodetemp;
		}
	}

	Int_t nmax = n+1;
	Int_t jmax = j+1;
	
	std::ifstream file;
	std::ofstream file2;
	file.open("RawConcentrations.txt");
	file2.open("Corrected Concentrations.txt",std::ofstream::app);
	
	Int_t run;
	Double_t concentration;
	Double_t counts;
	
	vector<Double_t> standard_counts(nmax);
	vector<Double_t> run_counts(jmax);
	
	vector<Double_t> run_concentration_averaged(jmax);
	vector<Double_t> run_concentration_corrected(run_files.size());
	vector<Double_t> run_concentration_corrected_index(run_files.size());
	
	vector<Double_t> standard_concentration_corrected(run_files.size());
	vector<Double_t> standard_concentration_corrected_index(run_files.size());
	vector<Double_t> standard_concentration_count(nmax);
	vector<Double_t> standard_concentration_averaged(nmax);

	while(1)
	{
		run = 0;
		counts=0;
		concentration=0;
		file >> run >> concentration >> counts;		
		
		for(n=0;n<nmax;n++)
		{
			Int_t index = 0;
			for(int f=0;f<n;f++)
			{
				index += standard_list[f].size();
			}
			for(int r=0;r<standard_list[n].size();r++)
			{
			if(run == standard_list[n][r])
			{
				Int_t ntot = standard_list[n].size();
				standard_concentration_corrected[index+r]=concentration/concentration_list[n];
				standard_concentration_corrected_index[index+r]=run;
				standard_concentration_averaged[n]+=concentration/ntot;
				standard_counts[n]+=counts;
				standard_concentration_count[n]=n;
			}
			}
		}
		
		for(j=0;j<jmax;j++)
		{
			Int_t index = 0;
			for(int f=0;f<j;f++)
			{
				index += run_list[f].size();
			}
			for(int r=0;r<run_list[j].size();r++)
			{
			if(run == run_list[j][r])
			{
				Int_t jtot = run_list[j].size();
				run_concentration_corrected[index+r] = concentration;
				run_concentration_corrected_index[index+r] = run;
				run_concentration_averaged[j]+=concentration/jtot;
				run_counts[j]+=counts;
			}
			}
		}
		
		if(!file.good()) break;
	}
	file.close();

	for(n=1;n<nmax;n++)
	{
		Int_t floor = *max_element(standard_list[n-1].begin(),standard_list[n-1].end());
		Int_t ceiling = *min_element(standard_list[n].begin(),standard_list[n].end());
		Int_t bottom = *min_element(standard_list[0].begin(),standard_list[0].end());
		Int_t top = *max_element(standard_list[nmax-1].begin(),standard_list[nmax-1].end());
		
		for(j=0;j<jmax;j++)
		{
			Int_t index = 0;
			for(int f=0;f<j;f++)
			{
				index += run_list[f].size();
			}
			Int_t getrunmax  = *max_element(run_list[j].begin(),run_list[j].end());
			
			if(getrunmax < bottom)
			{
				run_concentration_averaged[j] = run_concentration_averaged[j]*concentration_list[0]/standard_concentration_averaged[0];
				for(int r=0;r<run_list[j].size();r++)
				{
					run_concentration_corrected[index+r] = run_concentration_corrected[index+r]*concentration_list[0]/standard_concentration_averaged[0];
				}
			}
			else if(getrunmax > floor && getrunmax < ceiling)
			{
				run_concentration_averaged[j] = run_concentration_averaged[j]*0.5*((concentration_list[n-1]/standard_concentration_averaged[n-1])+(concentration_list[n]/standard_concentration_averaged[n]));
				for(int r=0;r<run_list[j].size();r++)
				{
					run_concentration_corrected[index+r] = run_concentration_corrected[index+r]*0.5*((concentration_list[n-1]/standard_concentration_averaged[n-1])+(concentration_list[n]/standard_concentration_averaged[n]));
				}
			}
			else if(getrunmax > top)
			{
				run_concentration_averaged[j] = run_concentration_averaged[j]*(concentration_list[nmax-1]/standard_concentration_averaged[nmax-1]);
				for(int r=0;r<run_list[j].size();r++)
				{
					run_concentration_corrected[index+r] = run_concentration_corrected[index+r]*(concentration_list[nmax-1]/standard_concentration_averaged[nmax-1]);
				}
			}
		}
	}
	
	Double_t final_concentration[100];
	Double_t final_error[100];
	Double_t run_number[100];
	Double_t run_error[100];
	
	for(j=0;j<jmax;j++)
	{
		if(run_counts[j]!=0)
		{
			Double_t error = (1/sqrt(run_counts[j]));
		}
		else Double_t error = 0;
	
		file2 << *min_element(run_list[j].begin(),run_list[j].end()) << "-" << *max_element(run_list[j].begin(),run_list[j].end()) << " " << run_name_list[j] << " " << run_concentration_averaged[j] << " " << error*run_concentration_averaged[j] << endl;
		concfinal[j]=run_concentration_averaged[j];
		final_error[j]=error*run_concentration_averaged[j];
		run_number[j]=j;
		run_error[j]=0;
	}
	file2.close();
	
	TCanvas* runs = new TCanvas("Runs","",1200,900);
	TGraphErrors* rungraph = new TGraphErrors(jmax,concfinal,runcount,errfinal,runerr1);
	rungraph->SetTitle("Sample Concentrations");
	rungraph->SetMarkerStyle(21);
	rungraph->GetXaxis()->SetLabelSize(.025);
	rungraph->GetYaxis()->SetTitle("Set Number");
	rungraph->GetYaxis()->CenterTitle();
	rungraph->GetYaxis()->SetTitleOffset(1.4);
	rungraph->GetYaxis()->SetLabelSize(.025);
	rungraph->SetMarkerSize(0.75);
	Double_t last = rungraph->GetXaxis()->GetXmax();
	rungraph->Draw("AP");
	auto legend = new TLegend(0.1,0.8,0.4,0.9);
	legend->AddEntry(rungraph,"Time-Averaged Concentrations","p");
	legend->Draw();
	 const char *text = "#bf{{}^{14}C /{}^{12}C Concentration}";
	TLatex latex = TLatex(0.5,0.02,text);
	latex.SetTextSize(0.035);
	Double_t size = latex.GetXsize();
	latex.SetTextAlign(11);
	latex.DrawLatexNDC(0.5-size/2,0.02,text);
	TText *t = new TText();
	t->SetTextSize(0.02);
	for(j=0;j<jmax;j++)
	{
	Double_t posx = concfinal[j];
	Double_t height = runcount[j]+0.2;
	t->DrawText(posx,height,run_name_list[j]);
	}
	runs->Update();

	TImage *img = TImage::Create();
	img->FromPad(runs);
	img->WriteImage("All Runs.png");
	
	TCanvas* runs2 = new TCanvas("Runs","",1200,900);
	TGraph* eachrun = new TGraph(totalruns,&concstandard_concentration_corrected_index[0],&concstandard_concentration_corrected[0]);
	runs2->SetGrid();
	eachrun->SetMarkerStyle(21);
	eachrun->GetXaxis()->SetNdivisions(run_files.size()/2);
	eachrun->GetXaxis()->SetRange(run_files[0],run_files[0]+run_files.size()-1);
	eachrun->GetXaxis()->SetLabelSize(.015);
	eachrun->Draw();
	
	TImage *img2 = TImage::Create();
	img2->FromPad(runs2);
	img2->WriteImage("Each Run.png");
	
	TCanvas* transmission = new TCanvas("Average Transmission","",1200,900);
	TGraph* transmission1 = new TGraph(nmax,&standard_concentration_count[0],&standard_concentration_averaged[0]);
	transmission1->SetMarkerStyle(21);
	transmission1->Draw();
	
	TImage *img3 = TImage::Create();
	img3->FromPad(transmission);
	img3->WriteImage("Average Transmission.png");
	
	TCanvas* transmission2 = new TCanvas("Run Transmission","",1200,900);
	TGraph* transmission3 = new TGraph(totalstans,&standard_concentration_corrected_index[0],&standard_concentration_corrected[0]);
	transmission3->SetMarkerStyle(21);
	transmission3->Draw();
	
	TImage *img4 = TImage::Create();
	img4->FromPad(transmission2);
	img4->WriteImage("Run Transmission.png");
}

void GetCathode(int number)
{
	name = cathode_name[number];
	conc = cathode_concentration[number];
}
