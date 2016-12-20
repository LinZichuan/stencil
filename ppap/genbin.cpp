#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sstream>
using namespace std;

int main() {
    string file = "case_720x360x38/A.txt";
    string binfile = "case_720x360x38/A.bin";
    FILE *fp;
    if ((fp = fopen(binfile.c_str(), "wb")) == NULL) {
        cout << "cannot open A.bin" << endl;
        return 0;
    }
    ifstream fin(file.c_str());
    int x,y,z,r;
    double v;
    timeval t1,t2;
    gettimeofday(&t1, NULL);
    int nx = 720, ny = 360, nz = 38;
    int size = nx*ny*nz*19;
    double *A = new double[size];
    int idx = 0;
    while (fin >> x >> y >> z >> r >> v) {
        A[idx] = v;
        idx++;
        //cout << x << endl;
        //cout << x << " " << y << " " << z << " " << r << " " << v << endl;
    }
    gettimeofday(&t2, NULL);
    float time = t2.tv_sec-t1.tv_sec;
    cout << "time = " << time << endl;
    fwrite(A, sizeof(double), size, fp);
    fclose(fp);
    fin.close();
    
    string fileb = "case_720x360x38/b.txt";
    string binfileb = "case_720x360x38/b.bin";
    FILE *fpb;
    if ((fpb = fopen(binfileb.c_str(), "wb")) == NULL) {
        cout << "cannot open b.bin" << endl;
        return 0;
    }
    ifstream finb(fileb.c_str());
    size = nx*ny*nz;
    double *b = new double[size];
    idx = 0;
    string strv;
    cout.precision(27);
    while (finb >> x >> y >> z >> strv) {
        stringstream ss;
        ss << strv;
        ss >> v;
        b[idx] = v;
        //cout << strv << endl;
        //cout << v << endl;
        idx++;
    }
    fwrite(b, sizeof(double), size, fpb);
    fclose(fpb);
    finb.close();

    string filex = "case_720x360x38/x0.txt";
    string binfilex = "case_720x360x38/x0.bin";
    FILE *fpx;
    if ((fpx = fopen(binfilex.c_str(), "wb")) == NULL) {
        cout << "cannot open x.bin" << endl;
        return 0;
    }
    ifstream finx(filex.c_str());
    size = nx*ny*nz;
    double *x0 = new double[size];
    idx = 0;
    while (finx >> x >> y >> z >> v) {
        x0[idx] = v;
        idx++;
    }
    fwrite(x0, sizeof(double), size, fpx);
    fclose(fpx);
    finx.close();

    return 0;
}
