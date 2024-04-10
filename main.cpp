#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <bitset>
#include <cfloat>


using namespace std;

ifstream fin("fisieralg.in");
ofstream fout("fisieralgout.out");


void incrucisare(vector<char> &cromozom1, vector<char> &cromozom2, int k){
    swap_ranges(cromozom1.begin()+k,cromozom1.end(), cromozom2.begin()+k);
}

bool mutatie(vector<char> &cromozom, double probabilitate){
    bool modificare = false;
    random_device rdm;
    mt19937  gen(rdm());
    uniform_real_distribution<double> dist4(0.0,1.0);
    for (int i_mut=0; i_mut<cromozom.size();i_mut++) {
        double probabilitate_aleatoare = dist4(gen);
        if (probabilitate_aleatoare < probabilitate) {
            modificare = true;
            if (cromozom[i_mut] == '0')
                cromozom[i_mut] = '1';
            else
                cromozom[i_mut] = '0';
        }

    }
    return modificare;
}


int cautareBinara(vector<double> &vec, double valoare ){
    int stanga = 0;
    int dreapta = vec.size() -1 ;

    while(stanga <= dreapta){
        int mijloc = stanga + (dreapta - stanga )/2;

        if(valoare >= vec[mijloc] && valoare < vec[mijloc+1])
            return mijloc+1;
        else if (valoare >= vec[mijloc+1] )
            stanga = mijloc+1;
        else
            dreapta = mijloc-1;
    }
    return 0;

}

int main() {
    int a, b, c;
    int interval_stanga;
    int interval_dreapta;
    int dimensiune_populatie;
    int precizie;
    double probabilitate_recombinare;
    double probabilitate_mutatie;
    int etape;
    double nr_random;
    double f_rec;
    double val_rec;
    int indiceInterval;

    fin >> a >> b >> c;
    fin >> interval_stanga;
    fin >> interval_dreapta;
    fin >> dimensiune_populatie;
    fin >> precizie;
    fin >> probabilitate_recombinare;
    fin >> probabilitate_mutatie;
    fin >> etape;

    int i;
    double discretizare;
    int biti = ceil(log2((interval_dreapta - interval_stanga) * pow(10, precizie)));
    discretizare = double((interval_dreapta - interval_stanga)) / pow(2.0, biti);

    //populatia initiala

    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<double> dist(interval_stanga, interval_dreapta);
    uniform_real_distribution<double> dist2(0, 1);
    uniform_int_distribution<> dist3(0, biti);

    vector<double> maxime;
    vector<double> medii;
    vector<double> functie;
    vector<double> x;
    vector<string> binar_x;
    vector<double> probabilitate_selectie;
    vector<double> functii_rec;
    vector<double> interval;
    vector<int> cromozomi_modificati;
    vector<int> ordine_selectare;
    vector<int> participanti;



    fout << "Populatia initiala" << endl;

    double suma = 0;

    for ( i = 1; i <= dimensiune_populatie; i++) {

                fout << i << ": ";
             nr_random = dist(gen);

             indiceInterval = ((nr_random - interval_stanga) / discretizare) + 1;

            bitset<22> bitiInterval(indiceInterval);

                fout << bitiInterval << " " << "x= " << nr_random << " f= " << fixed << setprecision(15)
                     << a * nr_random * nr_random + b * nr_random + c << endl;


            suma = suma + a * nr_random * nr_random + b * nr_random + c;

            functie.push_back(a * nr_random * nr_random + b * nr_random + c);
            x.push_back(nr_random);
            binar_x.push_back(bitiInterval.to_string());

        }

        //probabilitati selectie
    for(int etapa =1; etapa<=50; etapa++) {


        if(etapa==1)
            fout << endl << "Probabilitati selectie" << endl;
        for ( i = 0; i < dimensiune_populatie; i++) {
            if(etapa==1)
                fout << "cromozom  " << i + 1 << " probabilitate ";

            probabilitate_selectie.push_back(functie[i] / suma);
            if(etapa==1)
                fout << fixed << functie[i] / suma << endl;

        }

        if(etapa==1) {
            fout << endl;
            fout << "Intervale probabilitati selectie" << endl;
        }

        double initial = 0;
        double final = 0;
        interval.push_back(initial);
        for ( i = 0; i < dimensiune_populatie; i++) {
            if(etapa==1)
                fout << "[" << initial << " , ";
            final = final + probabilitate_selectie[i];
            initial = final;
            if(etapa==1)
                fout << final << " )" << endl;
            interval.push_back(final);

        }
        if(etapa==1)
            fout << endl;


        //generare aleatoare
        for ( i = 0; i < dimensiune_populatie; i++) {
            if(etapa==1)
                fout << "u= ";
             nr_random = dist2(gen);
            if(etapa==1)
                fout << nr_random << "  " << "selectam cromozomul " << cautareBinara(interval, nr_random) << endl;
            ordine_selectare.push_back(cautareBinara(interval, nr_random) - 1);
        }
        if(etapa==1)
            fout << endl;

        //dupa selectie
        if(etapa==1)
            fout << "Dupa selectie: " << endl;

        for ( i = 0; i < dimensiune_populatie; i++) {
            if(etapa==1) {
                fout << i + 1 << ": " << binar_x[ordine_selectare[i]] << " x= " << x[ordine_selectare[i]] << " f="
                     << functie[ordine_selectare[i]] << endl;
            }
        }
        if(etapa==1)
            fout << endl;

        //probabilitate incrucisare

        if(etapa==1)
            fout << "Probabilitatea de incrucisare " << probabilitate_recombinare << endl;

        for ( i = 0; i < dimensiune_populatie; i++) {
            nr_random = dist2(gen);
            if (nr_random < probabilitate_recombinare) {
                if(etapa==1) {
                    fout << i + 1 << ": " << binar_x[ordine_selectare[i]] << " u= " << nr_random << "<0.25 participa"
                         << endl;
                }
                participanti.push_back(i);
            } else
                if(etapa==1) {
                    fout << i + 1 << ": " << binar_x[ordine_selectare[i]] << " u= " << nr_random << endl;
                }
        }
        if(etapa==1)
            fout << endl;

        //daca participa nr impar ultimul nu mai este selectat

        if (participanti.size() % 2 == 1)
            participanti.pop_back();

        //recombinare si rezultat
        int j = 0;
        while (j < participanti.size()) {
            int punct_aleator = dist3(gen);
            if(etapa==1) {
                fout << "Recombinare dintre cromozomul " << participanti[j] + 1 << " cu cromozomul "
                     << participanti[++j] + 1
                     << ":" << endl;
            }
            vector<char> cromozom1(binar_x[ordine_selectare[participanti[j - 1]]].begin(),
                                   binar_x[ordine_selectare[participanti[j - 1]]].end());
            vector<char> cromozom2(binar_x[ordine_selectare[participanti[j]]].begin(),
                                   binar_x[ordine_selectare[participanti[j]]].end());
            if(etapa==1) {
                fout << binar_x[ordine_selectare[participanti[j - 1]]] << " "
                     << binar_x[ordine_selectare[participanti[j]]]
                     << " punct " << punct_aleator << endl;
            }

            incrucisare(cromozom1, cromozom2, punct_aleator);

            string cromozom1_modificat(cromozom1.begin(), cromozom1.end());
            string cromozom2_modificat(cromozom2.begin(), cromozom2.end());

            binar_x[ordine_selectare[participanti[j - 1]]] = cromozom1_modificat;
            binar_x[ordine_selectare[participanti[j]]] = cromozom2_modificat;
            j++;

            if(etapa==1)
                fout << " Rezultat  " << cromozom1_modificat << "  " << cromozom2_modificat << endl;


        }
        //dupa recombinare
        if(etapa==1) {
            fout << endl;
            fout << "Dupa recombinare: " << endl;
        }

        for ( i = 0; i < dimensiune_populatie; i++) {
            if(etapa==1)
                fout << i + 1 << ": " << binar_x[ordine_selectare[i]];
            int index = stoi(binar_x[ordine_selectare[i]], nullptr, 2);
             val_rec = interval_stanga + index * discretizare;
             f_rec = a * val_rec * val_rec + b * val_rec + c;

            if(etapa==1)
                fout << "  x =  " << val_rec << "  f= " << f_rec << endl;
        }

        //probabilitate mutatie
        if(etapa==1)
            fout << endl << "Probabilitate de mutatie pentru fiecare gena " << probabilitate_mutatie << endl;


        for ( i = 0; i < dimensiune_populatie; i++) {
            vector<char> cromozom_mut(binar_x[ordine_selectare[i]].begin(), binar_x[ordine_selectare[i]].end());
            if (mutatie(cromozom_mut, probabilitate_mutatie) == true) {
                cromozomi_modificati.push_back(i);
            }

            string cromozom_mut_modificat(cromozom_mut.begin(), cromozom_mut.end());

            binar_x[ordine_selectare[i]] = cromozom_mut_modificat;
        }
        //cromozomii care au fost modificati

        if (cromozomi_modificati.empty())
            if(etapa==1)
                fout << "Nu au existat modificari de cromozomi." << endl;

        else {
                if(etapa==1) {
                    fout << "Au fost modificati cromozomii: " << endl;
                    for (auto crom: cromozomi_modificati)
                        fout << crom + 1 << endl;
                }
        }

        //cromozomii dupa mutatie
        if(etapa==1)
            fout << " Dupa mutatie : " << endl;

        for ( i = 0; i < dimensiune_populatie; i++) {
            if(etapa==1)
                fout << i + 1 << ": " << binar_x[ordine_selectare[i]];

            int index = stoi(binar_x[ordine_selectare[i]], nullptr, 2);
            val_rec = interval_stanga + index * discretizare;
            f_rec = a * val_rec * val_rec + b * val_rec + c;

            if(etapa==1) {
                fout << "  x =  " << val_rec << "  f= " << f_rec << endl;

            }

            functii_rec.push_back(f_rec);

        }

        double maxim = DBL_MIN;

        for( auto functie_rec : functii_rec){
            if(functie_rec > maxim)
                maxim = functie_rec;
        }

        double suma_etapa1 = 0;
        for (auto valori_etapa1: functii_rec)
            suma_etapa1 = suma_etapa1 + valori_etapa1;
        double mean_etapa1 = suma_etapa1 / functii_rec.size();

        maxime.push_back(maxim);
        medii.push_back(mean_etapa1);

        functie.clear();
        x.clear();
        binar_x.clear();
        probabilitate_selectie.clear();
        interval.clear();
        ordine_selectare.clear();
        participanti.clear();
        cromozomi_modificati.clear();
        functii_rec.clear();


    }

  for( i =0; i <maxime.size(); i++){
      fout<<i+1<<": "<<"Maximul : "<<maxime[i]<<" si media : "<<medii[i]<<endl;
  }

    return 0;
}



