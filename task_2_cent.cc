//Invariant diff cross section vs pT for Pb-Pb Coll at 
//5.02TeV com energy for 0-5% Centrality
#include<iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TApplication.h"
#include "TFile.h"
#include "TLatex.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    TApplication theApp("hist", &argc, argv);
    Pythia pythia;

    pythia.readString("Beams:idA = 1000822080");
    pythia.readString("Beams:idB = 1000822080");
    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("HeavyIon:mode = 1");
    //pythia.readString("HeavyIon:Angantyr = on");
    pythia.readString("Beams:frameType = 1");

    // Initialize the Angantyr model to fit the total and semi-includive
    // cross sections in Pythia within some tolerance.
    pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.15");
    // These parameters are typicall suitable for sqrt(S_NN)=5TeV
    pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
    // A simple genetic algorithm is run for 20 generations to fit the
    // parameters.
    //double genlim[] = {2979.4, 2400.1, 1587.5, 1028.8, 669.9,
    //                 397.4, 220.3, 116.3, 54.5};
    // If you change any parameters these should also be changed.

    // The upper edge of the correponding percentiles:
   // double pclim[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

    pythia.readString("HeavyIon:SigFitNGen = 20");

    pythia.init();

    TH1F *pT_pion = new TH1F("pT_pion", "pT distribution of Pions", 1000, 0, 20.);
    TH1F *pT_kaon = new TH1F("pT_kaon", "pT distribution of Kaons", 1000, 0, 20.);
    TH1F *pT_proton = new TH1F("pT_proton", "pT distribution of Protons", 1000, 0, 20.);

    TH1F *nCharged = new TH1F("nCharged","Charged Multiplicity; Number of Charged Particles; Number of Events", 1000, 0, 50000.);

    int nEvents = 1000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        int multiplicity = 0;
        // Access the HIInfo object to get the centrality-related parameters
        HIInfo* hiInfo = pythia.info.hiInfo;
        double b = hiInfo->b(); // Impact parameter
        double Npart = hiInfo->nPartProj() +hiInfo->nPartTarg(); // Number of participants

        // Determine centrality based on Npart or b
        if (b <= 3.49 ) {// For Centrality 0 - 5%

            for (int i = 0; i < pythia.event.size(); ++i) {
                int id = pythia.event[i].id();
                double pT = pythia.event[i].pT();
                if (pythia.event[i].isCharged() && pythia.event[i].isFinal()){
                ++multiplicity;
            //if (pythia.event[i].isFinal()) {
                if ( id == 211 || id == -211) pT_pion->Fill(pT); //Pions
                if ( id == 321 || id == -321) pT_kaon->Fill(pT); // Kaons
                if (id == 2212 || id == -2212) pT_proton->Fill(pT); // Protons
                }
            }
            
            std::cout<<"Impact par  "<< b<<"   Npart   "<< Npart<<endl;
            std::cout<<"Multiplicity  "<< multiplicity<<endl;
        }
        nCharged->Fill(multiplicity); 
    }


    // Now convert the pT histogram to invariant yield: 1 / (2pi pT) * dN/dpT
    for (int i = 1; i <= pT_pion->GetNbinsX(); ++i) {
        double pT = pT_pion->GetBinCenter(i);  // Get the pT at the center of the bin

        // Avoid division by zero for pT = 0
        if (pT > 0) {
            double binContent_pion = pT_pion->GetBinContent(i);  // Get the bin content (dN/dpT)
            double Yield_pion = binContent_pion / (2 * TMath::Pi() * pT * nEvents);  // Invariant yield
            pT_pion->SetBinContent(i, Yield_pion);  // Set the bin content to the invariant yield
        }
    }

    for (int i = 1; i <= pT_kaon->GetNbinsX(); ++i) {
        double pT = pT_kaon->GetBinCenter(i);  // Get the pT at the center of the bin

        // Avoid division by zero for pT = 0
        if (pT > 0) {
            double binContent_kaon = pT_kaon->GetBinContent(i);  // Get the bin content (dN/dpT)
            double Yield_kaon = binContent_kaon / (2 * TMath::Pi() * pT * nEvents);  // Invariant yield
            pT_kaon->SetBinContent(i, Yield_kaon);  // Set the bin content to the invariant yield
        }
    }

    for (int i = 1; i <= pT_proton->GetNbinsX(); ++i) {
        double pT = pT_proton->GetBinCenter(i);  // Get the pT at the center of the bin

        // Avoid division by zero for pT = 0
        if (pT > 0) {
            double binContent_proton = pT_proton->GetBinContent(i);  // Get the bin content (dN/dpT)
            double Yield_proton = binContent_proton / (2 * TMath::Pi() * pT * nEvents);  // Invariant yield
            pT_proton->SetBinContent(i, Yield_proton);  // Set the bin content to the invariant yield
        }
    }
    TCanvas *c3 = new TCanvas("c3", "Charged Multiplicity in Pb-Pb Coll at 5.02 TeV", 800, 600);
    c3->SetLogy();
    nCharged->SetTitle("Charged Mutiplicity in Pb-Pb Collisions at 5.02 TeV");
    nCharged->Draw();

    TLatex *latex3 = new TLatex();
    latex3->SetNDC(); // Use normalized device coordinates (NDC) for positioning

    // Set text properties (optional)
    latex3->SetTextSize(0.03);
    latex3->SetTextAlign(23); // Align at top-left corner of the text box

    // Add the text annotations
    latex3->DrawLatex(0.5, 0.85, "#sqrt{s_{NN}} = 5.02 TeV"); 
    latex3->DrawLatex(0.5, 0.80, "Pb-Pb Collisions");
    latex3->DrawLatex(0.5, 0.75, "Impact Parameter b <= 3.49");
    c3->Update();

    TCanvas *c1 = new TCanvas("c1", "pT Distributions in Pb-Pb Collisions at 5.02 TeV", 800, 600);
    c1->SetLogy();
    pT_pion->SetMarkerColor(kRed);
    pT_kaon->SetMarkerColor(kBlue);
    pT_proton->SetMarkerColor(kGreen);
    
    pT_pion->SetTitle("pT Distribution in Pb-Pb Collisions at 5.02 TeV; pT (GeV); dN/dp_T;");
    pT_pion->SetMarkerStyle(20);
    pT_kaon->SetMarkerStyle(21);
    pT_proton->SetMarkerStyle(22);
    pT_pion->SetMarkerSize(0.75);
    pT_pion->Draw("AP");
    pT_kaon->Draw("P SAME");
    pT_proton->Draw("P SAME");

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(pT_pion, "Pions", "P");
    leg->AddEntry(pT_kaon, "Kaons", "P");
    leg->AddEntry(pT_proton, "Protons", "P");
    leg->Draw();

    // Add text annotations using TLatex
    TLatex *latex = new TLatex();
    latex->SetNDC(); // Use normalized device coordinates (NDC) for positioning

    // Set text properties (optional)
    latex->SetTextSize(0.03);
    latex->SetTextAlign(23); // Align at top-left corner of the text box

    // Add the text annotations
    latex->DrawLatex(0.6, 0.85, "#sqrt{s_{NN}} = 5.02 TeV");  // Center-of-mass energy
    latex->DrawLatex(0.6, 0.80, "Pb-Pb Collisions");         // Type of collision
    latex->DrawLatex(0.6, 0.75, Form("Events: %d", nEvents)); 
    latex->DrawLatex(0.6, 0.70, "Centrality: 0-5 %");

    c1->Update();
    TFile *outFile3 = new TFile("pT_distributions_Pb_cent.root", "RECREATE");
    c1->Write("mycanvas3");
    pT_pion->Write();
    pT_kaon->Write();
    pT_proton->Write(); 
    nCharged->Write();   
    outFile3->Close();

    pythia.stat();

    std::cout << "\nDouble-click on the canvas to quit." << std::endl;
    c1->WaitPrimitive();
    return 0;
}