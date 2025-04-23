#include <iostream>
#include <fstream>
#include <TMath.h>
#include <string>
#include <cmath>
#include <TF1.h>

using namespace std;

constexpr double pi = 3.14159265358979323846;
constexpr double c = 3e8; // m/s – prędkość światła


Double_t inter(Double_t *x, Double_t *p){
    double dx = x[0];

    // Parametry
    double A0        = p[0]; // amplituda głównego modu
    double A_plus1   = p[1]; // amplituda modu +1
    double A_minus1  = p[2]; // amplituda modu -1
    double Gamma     = p[3]; // szerokość spektralna
    double d         = p[4]; // długość rezonatora

    double A_plus1_tilde  = A_plus1 / A0;
    double A_minus1_tilde = A_minus1 / A0;

    double envelope = std::exp(-Gamma * std::abs(dx) / c);
    double cos1 = std::cos(pi * dx / d);
    double cos2 = std::cos(2 * pi * dx / d);

    double sqrt_term = 1
        + A_plus1_tilde * A_plus1_tilde
        + A_minus1_tilde * A_minus1_tilde
        + 2 * (A_plus1_tilde + A_minus1_tilde) * cos1
        + 2 * A_plus1_tilde * A_minus1_tilde * cos2;

    return (A0 / 3.0) * envelope * std::sqrt(sqrt_term);
}

Bool_t interfero(string filename = "data_corrected.txt"){
    TFile *output = new TFile("z27.root", "UPDATE");
    if(output->IsZombie()){return kFALSE;}

    TGraph *graph = new TGraph(filename.c_str(), "%lg %lg", "");
    Double_t fitMin = graph->GetXaxis()->GetXmin();
    Double_t fitMax = graph->GetXaxis()->GetXmax();

    TF1 *f = new TF1("interFit", inter, fitMin, 20, 5);
    f->SetParameters(1.15, 1.155, 1.395, 9e6, 15.7);

    f->SetParLimits(0, 0.3, 2.0);   // A0
    f->SetParLimits(1, 1, 2.0);   // A+1
    f->SetParLimits(2, 0.5, 2);   // A-1
    f->SetParLimits(3, 9e5, 2e9);
    f->SetParLimits(4, 10, 16);

    f->SetParName(0, "A0");        // amplituda głównego modu
    f->SetParName(1, "A_plus1");   // amplituda modu +1
    f->SetParName(2, "A_minus1");  // amplituda modu -1
    f->SetParName(3, "Gamma");     // szerokość spektralna
    f->SetParName(4, "d");         // długość rezonatora

    graph->Fit("interFit", "R");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    graph->SetMarkerStyle(8);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kOcean);
    graph->SetTitle("; #frac{1}{I_{0}} #left[#frac{cm^{2}}{W}#right]; #frac{1}{I_{1}} #left[#frac{cm^{2}}{W}#right]");
    graph->Draw("AP");
    graph->Write();

    return kTRUE;
}

Bool_t sim(){
    TFile *out = new TFile("output.root", "RECREATE");
    vector<Int_t> d = {10, 15, 20}; //cm
    Int_t G = 1.5e6; //Hz

    vector<Double_t> A_plus = {1, 1.2, 1.5};
    vector<Double_t> A_minus = {0.8, 1., 1.3};
    Double_t A_0 = 0.6;

    vector<TF1*> funcs;
    vector<TMultiGraph*> multi;
    vector<TLegend*> legends;
    
    
    for(int i = 0; i<3; i++){
        TMultiGraph *mg = new TMultiGraph(Form("resonatorLength%i", d[i]), Form("resonatorLength%i; #Delta x [cm]; V", d[i]));
        TLegend *leg = new TLegend();
        for(int j = 0; j<3; j++){
            TF1 *f = new TF1(Form("interFit_%i_%i", i, j), inter, -20, 70, 5);
            f->SetParameters(A_0, A_plus[j], A_minus[j], G, d[i]);
            funcs.push_back(f);

            TGraph *gr = new TGraph(f);
            gr->SetLineColor(92+j);
            leg->AddEntry(gr, Form("A(+1) = %.2f; A(-1) = %.2f", A_plus[j], A_minus[j]));
            mg->Add(gr);
        }
        legends.push_back(leg);
        multi.push_back(mg);
    }

    cout<<funcs.size()<<endl;

    // f->SetParName(0, "A0");        // amplituda głównego modu
    // f->SetParName(1, "A_plus1");   // amplituda modu +1
    // f->SetParName(2, "A_minus1");  // amplituda modu -1
    // f->SetParName(3, "Gamma");     // szerokość spektralna
    // f->SetParName(4, "d");         // długość rezonatora

    TCanvas *c1 = new TCanvas();
    c1->Divide(1, 3);
    for(int i = 0; i<3; i++){
        c1->cd(i+1);
        multi.at(i)->GetXaxis()->SetTitleSize(0.05);
        multi.at(i)->GetYaxis()->SetTitleSize(0.06);
        multi.at(i)->GetXaxis()->SetLabelSize(0.05);
        multi.at(i)->GetYaxis()->SetLabelSize(0.06);
        multi.at(i)->Draw("AL");
        legends.at(i)->Draw();
        multi.at(i)->Write();
    }

    c1->Write();
    
    return kTRUE;
}