//Invariant diff cross section vs pT for pp Coll at 5.02TeV CoM Energy
#include <iostream>
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TFile.h"
#include <TStopwatch.h>
#include "TLatex.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    TApplication theApp("hist", &argc, argv);
    TStopwatch timer;
    Pythia pythia;
    
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("HardQCD:all = on");
    pythia.readString("SoftQCD:all = on"); 
    pythia.init();

    TH1F *pT_pion = new TH1F("pT_pion", "pT distribution of Pions", 1000, 0, 20);
    TH1F *pT_kaon = new TH1F("pT_kaon", "pT distribution of Kaons", 1000, 0, 20);
    TH1F *pT_proton = new TH1F("pT_proton", "pT distribution of Protons", 1000, 0, 20);
    
    TH1F *nCharged = new TH1F("nCharged","Charged Multiplicity; Number of Charged Particles; Number of Events", 1000, 0, 500.);

    timer.Start();
    // Generate events and fill histograms
    int nEvents = 10000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        int multiplicity = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
                int id = pythia.event[i].id();
                double pT = pythia.event[i].pT();
                double eta = pythia.event[i].eta();
                if (pythia.event[i].isCharged() && pythia.event[i].isFinal() && eta >=-2.0 && eta<=2.0){
                ++multiplicity;
                if ( id == 211 || id == -211) pT_pion->Fill(pT); //Pions
                if ( id == 321 || id == -321) pT_kaon->Fill(pT); // Kaons
                if (id == 2212 || id == -2212) pT_proton->Fill(pT); // Protons
                }
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
    timer.Stop();

    TCanvas *c3 = new TCanvas("c3", "Charged Multiplicity in pp Coll at 5.02 TeV", 800, 600);
    c3->SetLogy();
    nCharged->SetTitle("Charged Mutiplicity in pp Collisions at 5.02 TeV");
    nCharged->Draw();

    TCanvas *c1 = new TCanvas("c1", "pT Distributions", 800, 600);
    c1->SetLogy();
    pT_pion->SetMarkerColor(kRed);
    pT_kaon->SetMarkerColor(kBlue);
    pT_proton->SetMarkerColor(kGreen);
    pT_pion->SetTitle("pT Distribution; p_{T} (GeV); #frac{1}{2#pi p_{T}} #frac{dN}{dp_{T}};");
    pT_pion->SetMarkerStyle(20);
    pT_kaon->SetMarkerStyle(21);
    pT_proton->SetMarkerStyle(22);
    pT_pion->SetMarkerSize(0.62);
    pT_pion->Draw("P");
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
    latex->DrawLatex(0.6, 0.85, "#sqrt{s} = 5.02 TeV");  // Center-of-mass energy
    latex->DrawLatex(0.6, 0.80, "pp Collisions");         // Type of collision
    latex->DrawLatex(0.6, 0.75, Form("Events: %d", nEvents)); 

    c1->Update();
    // Save the histograms to a file
    TFile *outFile1 = new TFile("pT_distributions.root", "RECREATE");
    c1->Write("mycanvas1");
    pT_pion->Write();
    pT_kaon->Write();
    pT_proton->Write();
    nCharged->Write(); 
    outFile1->Close();

    //c1->SaveAs("InvariantDiffCrossSection_pions.png");
    //std::cout << "Real time: " << timer.RealTime() << " seconds" << std::endl;
    //std::cout << "CPU time: " << timer.CpuTime() << " seconds" << std::endl;
    timer.Print();
    std::cout << "\nDouble-click on the canvas to quit." << std::endl;
    c1->WaitPrimitive();
    return 0;
}
