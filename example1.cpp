#include <iostream>
#include <fstream>
#include <sstream>
#include "lrmodel.h"
#include "lrfaxial.h"
#include "compress.h"
#include "lrfio.h"
#include "bspline123d.h"
#include "bsfit123.h"
#include "reconstructor.h"
#include <cmath>

int main()
{
// 1. Create a square 8x8 sensor array with 4.21 mm pitch
    LRModel lrm(64);
    double step = 4.21;
    double shift = step*3.5;
    for (int i=0; i<64; i++) {
        double y = -(i/8 * step - shift);
        double x = i%8 * step - shift;
        lrm.AddSensor(i, x, y);
//        std::cout << x << ", " << y << std::endl;
    }

// 2. Group sensors assuming square symmetry
    lrm.MakeGroupsSquare();

    std::cout << "Number of Sensors: " << lrm.GetSensorCount() << std::endl;
    std::cout << "Number of Groups: " << lrm.GetGroupCount() << std::endl << std::endl;

    for (int i=0; i<lrm.GetGroupCount(); i++) {
        std::cout << "Sensors in Group #" << i << ": ";
        for (auto q : lrm.GroupMembers(i))
            std::cout << q << ", ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

// 3. Create default LRF: axial with 10 spline nodes
// Rmax doesn't matter at this point, it will be auto-adjusted before fitting
    LRFaxial *mylrf = new LRFaxial(42., 10);
// Compression: dual slope with k=10, r0=7 and lambda=4
    DualSlopeCompress *compr = new DualSlopeCompress(10., 7., 4.);
    mylrf->SetCompression(compr);
// Add constraints so that LRFs always get nice shape
    mylrf->setNonNegative(true);
    mylrf->SetNonIncreasing(true);
    mylrf->SetFlatTop(true);

    for (int i=0; i<lrm.GetGroupCount(); i++) {
        lrm.SetGroupLRF(i, mylrf->clone());
        ((LRFaxial*)(lrm.GetGroupLRF(i)))->SetOrigin(lrm.GetGroupX(i), lrm.GetGroupY(i));
    }

// 4. Load simulated flood
    // format: a0, ... a63, nPhotons, x, y
    std::vector <std::vector <double> > Data;

    std::ifstream f("Simulation_10k.txt");
    if (!f.good()) {
        std::cout << "Unpack data/Simulation_10k.txt.zip into work directory first" << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(f, line)) {
        std::vector <double> evt;
        std::istringstream iss(line); //put line into stringstream
        double val;
        while(iss >> val) //read word by word
            evt.push_back(val);
        Data.push_back(evt);
        evt.clear();
    }

// 5. Fit the LRFs to the flood data
    std::vector < LRFdata > d0;
    for (auto q : Data) {
        d0.push_back(LRFdata({q[65], q[66], 0, q[0]}));
    }

// adjust Rmax
    for (int i=0; i<lrm.GetGroupCount(); i++) {
        LRFaxial *lrf = dynamic_cast<LRFaxial*> (lrm.GetGroupLRF(i));
        lrf->SetRmax(lrm.GetGroupMaxR(i, d0));
//        std::cout << lrm.GetGroupMaxR(i, d0) << std::endl;
    }

    for (int j=0; j<64; j++) {
        for (int i=0; i<d0.size(); i++) {
            d0[i][3] = Data[i][j];
        }
        lrm.AddFitData(j, d0);
    }

    for (int i=0; i<lrm.GetGroupCount(); i++)
        lrm.FitGroup(i);

    std::cout << std::endl << "--------------------------------------------------" << std::endl << std::endl;

// 6. Reconstruct the flood events

    Reconstructor reco(&lrm);
    reco.InitMinimizer();
    reco.setCogRelCutoff(0.1);
    reco.setEnergyCalibration(0.005);

    std::vector <bool> sat(64, false); // no saturation

// begin by reconstructing and printing results for the first 10 events
    std::cout << "Here are the first 10 recnstructed events:" << std::endl;
    for (int i=0; i<10; i++) {
        reco.ProcessEvent(Data[i], sat);
        std::cout << "Event " << i << " status=" << reco.getRecStatus() << std::endl;
        std::cout << "X " << Data[i][65] << " vs " << reco.getRecX() << std::endl;
        std::cout << "Y " << Data[i][66] << " vs " << reco.getRecY() << std::endl;
        std::cout << "E " << reco.getRecE() << std::endl;
        std::cout << "Eguess " << reco.getGuessE() << std::endl << std::endl;
    }

// and then reconstruct the whole flood and dump the result into a file
    std::ofstream recfile;
    recfile.open ("reconstruction.txt", std::ios_base::app);
    for (auto d : Data) {
        reco.ProcessEvent(d, sat);
        recfile << d[65] << " " << reco.getRecX() << " " << d[66] << " " << reco.getRecY() << std::endl;
    }
    recfile.close();

// 7. Save the detector model as JSON
    std::ofstream datafile;
    datafile.open ("LRM_square8x8.json", std::ios_base::app);
    datafile << lrm.GetJsonString();
    datafile.close();

    return 0;
}
