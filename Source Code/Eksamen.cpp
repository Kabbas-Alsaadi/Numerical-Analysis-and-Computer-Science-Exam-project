#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <random>
#include <stdlib.h>

#define nmax 100
#define eps 0.000000000000000000000000001
#define MTIT 100

using namespace std;


void Brugerorientering();

int vaelg_A(double A[nmax][nmax]);
void hent_matrix(double M[nmax][nmax],int &n, string fn);
void udskriv_matrix(double A[nmax][nmax], int n, int m);
void udskriv_Totalmatrix(double A[nmax][nmax+1], int n, int m);
void udskriv_vektor(double v[nmax], int n);
void gem_Aogb(double A[nmax][nmax], double b_vektor[nmax],int n, int m, string fn);
void hent_Aogb(double A[nmax][nmax], double b[nmax],int &n, int &m, string fn);

void Valgb(double A[nmax][nmax], double b[nmax], int &n, int &m);


////2) x* findes ved modificeret Gram-Schmidt
void MinKvad_ved_Mod_GS(double A_matrix[nmax][nmax], double b_vektor[nmax], int n, int m);
void Kopivektorimatrix(double vektor[nmax],int n, int j, double Resultatmatrix[nmax][nmax]);
double Prikprodukt(double qNy[nmax], double Agl[nmax][nmax], int n, int k);
double Søjlelængde_i_matrix(double Matrix[nmax][nmax], int n, int k);
double Vektorlængde(double vektor[nmax], int n);


////Frembringelse af datagrundlag C
void Inddelingafr(double r[nmax],double xk[7], int n);
void NewtonRaphsonFordelingskontrol(double xk[7]);
double FixByTaylor(double x);
double fmrk(double x);
void StakSøjler(double vektor[nmax], double Matrix[nmax][nmax], int n, int m);
void Ydreprodukt_af_vektorer(double vektor_a[nmax], double vektor_b[nmax], double Resultatmatrix[nmax][nmax], int n, int m);
void Simuler_r(double r[nmax], int n);
void Datagrundlag_c_1(double x[nmax],double &deltax,int &mx,double emx[nmax],double y[nmax],double &deltay,int &ny,double eny[nmax],double r[nmax], int &n, int &m);
void Dan_A_og_b(double x[nmax],double deltax,int mx,double emx[nmax],double y[nmax],double deltay,int ny,double eny[nmax],double r[nmax], int n, int m, double A_matrix[nmax][nmax], double b_vektor[nmax]);

////Generelle hjælpefunktioner
void Transponer(double Matrix[nmax][nmax],double Transponeret[nmax][nmax], int n, int m);
void Matrixvektorprodukt(double Matrix[nmax][nmax], double vektor[nmax], double Resultatvektor[nmax], int n);
void MatrixMatrixprodukt(double Matrix1[nmax][nmax], double Matrix2[nmax][nmax], double Resultatmatrix[nmax][nmax], int n, int m);
void Matrixvektorprodukt_med_m(double Matrix[nmax][nmax], double vektor[nmax], double Resultatvektor[nmax], int n, int m);

////Frembringelse datagrundlag d
void Hent_A_til_d(double A[nmax][nmax], int n, int m, string fn);
void Datagrundlag_d(double A_matrix[nmax][nmax], double b_vektor[nmax], int &n, int &m);


////Rækkereduktion
void Backsub_uden_totalmatrix(double A_matrix[nmax][nmax], double b_vektor[nmax], int n, double Løsning[nmax]);
void Backsub(double Totalma[nmax][nmax], double x[nmax], int n);
void DelvisPivotering(double ma[nmax][nmax], int j,  int n);
void Gauss(double ma[nmax][nmax], int n);
void LineærLøsning(double Totalma[nmax][nmax], double Løsning[nmax], int n);
void LineærLøsning_uden_totalmatrix(double A_matrix[nmax][nmax], double b_vektor[nmax], int n, double Løsning[nmax]);


////4): Find egenløsninger
void Test_Egenløsning();
void Egenløsninger(double ATA_matrix[nmax][nmax], int m, double Lambda[nmax], double Egenvektorer[nmax][nmax]);
double LgdDifVek(double v1[nmax], double v2[nmax], int n);
void NormerVektor(double v[nmax], double vnorm[nmax], double &L, int n);
void Kopicib(double c[nmax], int n, double b[nmax]);
void Pk_kalkulation(double Pgl[nmax][nmax], double Egenvektorer[nmax][nmax], double Lny, int k, int m, double Pny[nmax][nmax]);

////3) Bestem x* ved min f(x)
void MinKvad_ved_minf(double A_matrix[nmax][nmax],double ATA_matrix[nmax][nmax],double b_vektor[nmax], double ATb_vektor[nmax], double Lambda[nmax], double Egenvektorer[nmax][nmax], int n, int m);
double T_InkredsningsAlgoritme(double A_matrix[nmax][nmax], double b_vektor[nmax], double x[nmax], int n, int m, double D, double Epsilon_GS, double sk[nmax]);
double f(double A_matrix[nmax][nmax], double b_vektor[nmax], double x[nmax], int n, int m);
void Fremstil_Plot(double A_matrix[nmax][nmax], double b_vektor[nmax], double x0[nmax], double s0[nmax], int n, int m);
int main() {
	char Valgabcd, ValgGemB, ValgGemD, ValgKontrolr, Valg124, Valg3, ValgIgen, ValgEgen;
	double A_matrix[nmax][nmax], b_vektor[nmax], AT_matrix[nmax][nmax], ATA_matrix[nmax][nmax], ATb_vektor[nmax];
	double MinKvadx[nmax], Lambda[nmax], Egenvektorer[nmax][nmax];
	int n, m;
	double xk_til_Fordelingskontrol[nmax];

	////Data grundlag c variable
	double x[nmax], deltax, emx[nmax],y[nmax],deltay,eny[nmax],r[nmax];
	int ny, mx;

	string a__IndlæsAogb, gemAogb__b;
do
{
	Brugerorientering();
	cout << "Vælg a, b, c eller d her: "; cin >> Valgabcd;

	switch(Valgabcd)
	{
	case 'a':
			{
			cout << "Du indtastede a, dermed valgte du a): Indlæs A, n, m, og vektoren b fra en prækonstrueret datafil " << endl;
			cout << "På dette sted i programmet udføres: Indlæsning af A, n, m og b fra en prækonstrueret datafil." << endl;
			hent_Aogb(A_matrix, b_vektor, n, m, "IndlæsAogB.txt");
			cout << "Matricen A ser ud således: " <<endl << endl;
			udskriv_matrix(A_matrix, n, m);
			cout << endl << endl << "vektoren b ser ud således: "<< endl << endl;
			udskriv_vektor(b_vektor, n);
			cout<< "Test n:"<<endl;
			}break;

	case 'b':
			{
				cout << "Du indtastede b - dermed valgte du b): Du indtaster selv A, n, m og b" << endl;
				cout << "På dette sted i programmet indtaster du selv A, n, m og b" << endl;
				Valgb(A_matrix, b_vektor, n, m);
				cout << "Du har nu indtastet n, m, elementer for A og elementerne for b. Ønsker du at gemme dette?" << endl;
				cout <<"Indtast j her for at gemme: "; cin >> ValgGemB;
				if(ValgGemB == 'j')
				{
					gem_Aogb(A_matrix, b_vektor, n, m, "Abb.txt");
				cout << "Programmet har nu gemt dine indtastede værdier for n og m, og for samtlige elementer i A og b" << endl;

				}
				else
				{
				cout << "Du har valgt ikke at gemme dine indtastede værdier for n og m, og for samtlige elementer i A og b. Vi går videre" << endl;
				}

			}break;


	case 'c':
			{
			cout << "Du indtastede c, dermed valgte du c): ";
			cout << "Simulér målinger (x_i, y_i, z_i), hvor i = 1, 2,...,n, for en multipel regression:" << endl;
			cout << "z_i = a*x_i + b*y_i + c + eps_i, hvor eps_i er simulerede observationer af N(0, sigma)^2 og deraf dannes A og b" << endl;
			cout << "På dette sted i programmet simuleres målinger for en multipel regression som skal bruges til at danne A og b" << endl << endl;
			Datagrundlag_c_1(x, deltax, mx, emx,y,deltay,ny,eny,r,  n, m);
			cout << "Du har nu muligheden for at programmet skal udføre en fordelingskontrol. Ønsker du at programmet skal gøre dette, indtastes j her: " << endl; cin >> ValgKontrolr;
				if(ValgKontrolr == 'j')
				{
				cout << "Du har valgt at programmet skal udføre en fordelingskontrol: " << endl;
				NewtonRaphsonFordelingskontrol(xk_til_Fordelingskontrol);
				Inddelingafr(r,xk_til_Fordelingskontrol, n);
				cout << "Programmet har nu udført en fordelingskontrol, vi går videre" << endl;
				}
				else
				{
				cout << "Du har valgt at programmet ikke skal udføre en fordelingskontrol, vi går videre" << endl;
				}
			Dan_A_og_b(x, deltax, mx, emx, y, deltay, ny, eny, r, n, m, A_matrix, b_vektor);

 			}break;

	case 'd':
			{
			cout << "Du indtastede d, dermed valgte du d): Dan A og b så A*x=b har en brugervalgt mindste kvadraters løsning x*" << endl;
			cout << "På dette sted i programmet vælger du x*, hvorefter programmet danner A og b således A*x=b er konsistent" << endl;
			Datagrundlag_d(A_matrix, b_vektor, n, m);
			cout << "Ønsker du at gemme A, n, m og b? Hvis ja, så indtast j" << endl; cin >> ValgGemD;
				if(ValgGemD == 'j')
				{
					cout << "Du har valgt at gemme A, n, m og b. Vi går videre" << endl;
					gem_Aogb(A_matrix, b_vektor, n, m, "Abd.txt");
				}
				else
				{
					cout << "Du har valgt ikke at gemme A, n, m og b. vi går videre" << endl;
				}
			}break;
	}

	Transponer(A_matrix, AT_matrix, n, m);

	MatrixMatrixprodukt(AT_matrix, A_matrix, ATA_matrix,n, m);

	cout << endl << endl <<"A matrix: " << endl;

	udskriv_matrix(A_matrix, n, m);

	cout<<endl<<endl<<"AT matrix: "<<endl;
	udskriv_matrix(AT_matrix, m, n);

	cout << endl<<endl << "ATA matrix: "<<endl;

	udskriv_matrix(ATA_matrix, m, m);

	Matrixvektorprodukt(AT_matrix, b_vektor, ATb_vektor, n);

	cout << endl << endl << "ATb vektor: " << endl;

	udskriv_vektor(ATb_vektor, m);

	cout << endl << endl;
	cout << "Du har fået programmet til at frembringe datag"
			"rundlaget ved valgmulighed: " << Valgabcd << ". Du kan nu vælge hvad programmet skal" << endl;
	cout << "gøre med dette datagrundlag. Der kan ske en bestemmelse af mindste kvadraters løsning x*, " << endl;
	cout << "til et overbestemt system af lineære ligninger A*x=b, eller en bestemmelse af egenløsninger " << endl;
	cout << "(lambda_k, v_k), hvor k = 1, 2, ...,n, til en symmetrisk matrix A. " << endl << endl;
	cout << "Bestemmelse af mindste kvadraters løsning x* kan du vælge at skulle foretages ved 3 forskellige metoder:" << endl << endl;
	cout << "1) x* findes som løsning til normalligningerne: Transpose(A)*A*x=Transpose(A)*b" << endl;
	cout << "2) x* findes ved at udfører først QR-faktorisering af A ved modificeret Gram-Schmidt ortogonalisering og dernæst udføres." << endl;
	cout << "løses R*x = Transpose(Q)*b ved backwards substitution" << endl;
	cout << "Dernæst findes x* som løsning til R*x=Transpose(Q)*b" << endl;
	cout << "3) x* findes ved minimering af f(x) = (Det(b-A*x))^2. BEMÆRK: Denn kan du kun vælge hvis du vælger 4 først" << endl << endl;
	cout << "Bestemmels af egenløsningerne kan findes ved én metode:" << endl;
	cout << "4) (lambda_k, v_k) findes ved en udvidet version af potensmetoden" << endl << endl;
	cout << "Ønskes 1, da indtastes 1, ønskes 2, da indtastes 2, ønskes 4, da indtast 4. Indtast venligst her: "; cin >>Valg124;


	switch(Valg124)
	{
	case '1':
			{
			cout << "Du indtastede 1 - dermed valgte du 1): x* findes som løsning til normalligningerne: Transpose(A)*A*x=Transpose(A)*b" << endl;
			cout << "På dette sted i programmet findes den mindste kvadraters løsning som en løsning til normalligninger: Transpose(A)*A*x=Transpose(A)*b " << endl;


			LineærLøsning_uden_totalmatrix(ATA_matrix, ATb_vektor, m, MinKvadx);


			cout << endl << endl<< "Mindste kvadraters løsningen x* er lig med: " << endl;
			udskriv_vektor(MinKvadx, m);
			}break;


	case '2':
			{
			cout << "Du indtastede 2 - dermed valgte du 2): x* findes ved at udfører først QR-faktorisering af A ved modificeret Gram-Schmidt ortogonalisering og dernæst udføres." << endl;
			cout << "løses R*x = Transpose(Q)*b ved backwards substitution" << endl;
			cout << "På dette sted i programmet findes x* ved at udføre backwards substitution på en QR-faktorisering A-matrix med modificeret Gram-Schmidt ortogonalisering" << endl;
			cout <<"Her ses A matricen: "<<endl; udskriv_matrix(A_matrix,n,m);
			MinKvad_ved_Mod_GS(A_matrix, b_vektor, n, m);
			}break;

	case '4':
			{
			cout << "Du indtastede 4 - dermed valgte du 4): (lambda_k, v_k) findes ved en udvidet version af potensmetoden" << endl;
			cout << "På dette sted i programmet findes egenløsninger matricen Transpose(A)*A ved en udvidet version af potensmetoden" << endl<< endl;
			cout << "((Programmet finder egenløsningerne til matricen Transpose(A)*A))" << endl << endl;
			Egenløsninger(ATA_matrix, m, Lambda, Egenvektorer);
			cout << endl << endl << "Her ses egenværdierne: " << endl; udskriv_vektor(Lambda, m);
			cout << endl << endl << "Her ses Egenvektorerne: " << endl; udskriv_matrix(Egenvektorer, m, m);

			cout << "Ønsker du at teste del af programmet i forhold til at finde egenløsninger for en matrix med kendte egenløsninger?"<<endl;cin>>ValgEgen;

			if(ValgEgen == 'j')
			{
				Test_Egenløsning();
			}
			cout << "Du kan nu vælge om du vil finde mindste kvadraters løsning x* ved minimering af funktionen f(x) = Transpose(b-A*x)*(b-A*x)" << endl;
			cout << "For at finde x* ved denne metode, indtast da j her: "; cin >> Valg3;

			if(Valg3 == 'j')
			{
				cout << "Du har valgt at programmet skal finde x* ved at minimere funktionen f(x) = Transpose(b-A*x)*(b-A*x)"<< endl;
				MinKvad_ved_minf(A_matrix,ATA_matrix,b_vektor, ATb_vektor, Lambda, Egenvektorer,n, m);
			}
			else
			{
				cout << "Du har valgt at programmer ikke skal finde x* ved at minimere funktionen f(x). Vi går videre" << endl;
			}
			}break;

	}
	cout << "Du er nu nået til slutningen af programmet. Ønsker du at køre programmet igen. Hvis det ønskes, så indtast j her:"; cin >> ValgIgen;
}
while(ValgIgen == 'j');

cout << "Tak for du har benyttet dette program." "Over and out.";
	return 0;
}

void Test_Egenløsning()
{
	double Test_A[nmax][nmax], D[nmax], Test_Lambda[nmax], S[nmax][nmax], ST[nmax][nmax], M[nmax][nmax], S_Lambda[nmax][nmax];
	double Egenvektorer[nmax][nmax], Lambda[nmax];
	int n = vaelg_A(Test_A);
	cout<<endl<<endl<<"Test_A: "<<endl; udskriv_matrix(Test_A, n, n);
	for(int k=0;k<n;k++)D[k]=Søjlelængde_i_matrix(Test_A, n, k);
	cout<<endl<<endl<<"Længden på søjlerne af Test_A: "<<endl; udskriv_vektor(D, n);
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			S[j][i]= Test_A[j][i]/D[i];
		}
	}
	Transponer(S, ST, n, n);
	cout<< "Her ses testmatricens egenvektorer:"<<endl; udskriv_matrix(S, n, n);
	cout << "Du skal nu indtaste egenværdierne for matricen." <<endl<<endl;
	for(int i=0;i<n;i++)
	{
		cout<< "Indtast venligst den "<<i+1<<". egenværdi: ";cin>>Test_Lambda[i]; cout<<endl;
	}

	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			S_Lambda[j][i]= S[j][i]*Test_Lambda[i];
		}
	}
	MatrixMatrixprodukt(S_Lambda, ST, M, n, n);
	cout<<"Matricen med de kendte egenløsningener ser ud således: "<<endl; udskriv_matrix(M, n, n);
	cout<<"Nu findes egenløsninger ved Pn kalkulation."<<endl;

	Egenløsninger(M, n , Lambda, Egenvektorer);
	cout<<"Her ses egenværdierne som blev brugt til at konstruere testmatricen: "<<endl;udskriv_vektor(Test_Lambda, n);
	cout<<"Her ses testmatricen egenvektorer: "<<endl; udskriv_matrix(S, n, n);
	cout<<endl<<endl<<"Her ses egenværdierne ved Pn kalkulation: "<<endl; udskriv_vektor(Lambda, n);
	cout<<"Her ses egenvektorene ved Pn kalkulation: "<<endl; udskriv_matrix(Egenvektorer, n, n);

}
////3) Bestemmelse af x* ved min f(x), kan være fejl med valg af m og/rllrt n, skiftet rundt eller forkert
void MinKvad_ved_minf(double A_matrix[nmax][nmax],double ATA_matrix[nmax][nmax],double b_vektor[nmax], double ATb_vektor[nmax], double Lambda[nmax], double Egenvektorer[nmax][nmax], int n, int m)
{
	double I_Hesse_matrix[nmax][nmax], I_H_lambda[nmax], PD[nmax][nmax], PT[nmax][nmax],gradf[nmax],D, Epsilon, Epsilon_GS,N;
	double ATAx[nmax],s0[nmax],sk[nmax],x0[nmax], xGl[nmax], xNy[nmax], vk[nmax], LgdDiff,gradf_Lgd; LgdDiff=0;gradf_Lgd=0;

	double t;t=0;
	int k=0;
	for(int i=0; i<m; i++) I_H_lambda[i]= 1/(2*Lambda[i]);
	for(int i=0; i<m; i++){for(int j=0;j<m;j++){PD[j][i]=Egenvektorer[j][i]*I_H_lambda[i];}}
	Transponer(Egenvektorer, PT, m, m);
	MatrixMatrixprodukt(PD, PT, I_Hesse_matrix, m, m);
	cout<<endl<<endl<<"Her ses den inverse Hesse-matrix for Transpose(A)*A: "<<endl;udskriv_matrix(I_Hesse_matrix, m, m);
	cout << "Du skal nu indtaste værdier for elementerne i x0, som er startvektoren for iterationen."<<endl;
	for(int i=0;i<m;i++){cout<<"Indtast venligst det "<<i+1<<". element: ";cin>>x0[i];cout<<endl; xNy[i]=x0[i]; xGl[i]=x0[i];}
	cout <<endl<<"Indtast venligst tolerancen for stopkriteriet i algoritmen Epsilon_1: ";cin>>Epsilon;
	cout<<endl<<"Indtast venligst det maksimal tilladte iterationer for algoritmen, N: ";cin>>N;
	cout<<endl<<"Indtast venligst tolerancen for stopkriteriet i under-algoritmen 'Golden section': ";cin>>Epsilon_GS;
	cout<<endl<<"Indtast venligst startsteplængden for indkredsningsalgoritmen, som finder et interval for Golden section: ";cin>>D;
	for(int i=0;i<m;i++){for(int j=0; j<m;j++){I_Hesse_matrix[i][j]=I_Hesse_matrix[i][j]*(-1);}};

	cout<<endl<<endl<<"her ses funktionsværdien for f(x0) = "<<f(A_matrix, b_vektor, x0, n,m);
	for(int i=0; i<m;i++) ATb_vektor[i]=ATb_vektor[i]*2;
	do
	{
		//Finder gradf

		for(int i=0;i<m;i++)ATAx[i]=0;
		Matrixvektorprodukt(ATA_matrix, xGl,ATAx, m);

		for(int i=0;i<m;i++) ATAx[i]=ATAx[i]*2;
		for(int i=0;i<m;i++) gradf[i]=ATAx[i]-ATb_vektor[i];

		Matrixvektorprodukt(I_Hesse_matrix, gradf, vk, m);

		for(int i=0; i<m;i++) sk[i]=vk[i]/Vektorlængde(vk, m);if(k==0){ for(int i=0; i<m; i++){ s0[i]=sk[i];}}



		t=T_InkredsningsAlgoritme(A_matrix, b_vektor, xGl, n, m, D, Epsilon_GS, sk);


		for(int i=0; i<m;i++)xNy[i]=xGl[i]+t*sk[i];

		t=0;

		LgdDiff=LgdDifVek(xNy, xGl, m);
		gradf_Lgd=Vektorlængde(gradf, m);
		for(int i=0;i<m;i++)xGl[i]=xNy[i];
		k++;
	}while(LgdDiff>Epsilon && gradf_Lgd>Epsilon && k<N);
	cout<<endl<<endl<<"Her ses mindste kvadraters løsningen x*: "<<endl;udskriv_vektor(xNy, m);
	Fremstil_Plot(A_matrix, b_vektor,x0, s0, n, m);

}

void Fremstil_Plot(double A_matrix[nmax][nmax], double b_vektor[nmax], double x0[nmax], double s0[nmax], int n, int m)
{
	int k=1000;
	double t[k+1], funktionsværdi[k+1], x_t_sk[nmax][k+1], Temp[nmax], t_L;
	string fn = "Plot.txt";

	cout<<"Denne del af programmet skaber et datagrundlag for et plot som direkte kan kopieres over i"<<endl;
	cout<<"i Mathematica for at plotte f(x) som funktion af t: (t, (f(x0(t)))), så man kan aflæse minimum derfra"<<endl;
	cout<<"Der plottes i intervallet: 0 <= t <= t_L"<<endl;

	cout<<endl<<"Du bedes indtaste t_L, som er den største værdi t kan antage i plottet: "<<endl;cin>>t_L;


	t[k]=t_L;
	t[0]=0;


	for(int i=1;i<k+1;i++) t[i]=(t_L/k)*i;

	for(int i=0;i<k+1;i++)
	{
		for(int j=0;j<n;j++)
		{
			x_t_sk[j][i]=x0[j]+t[i]*s0[j];
		}
	}

	for(int i=0;i<k+1;i++)
	{
		for(int j=0;j<n;j++)
		{
			Temp[j]=x_t_sk[j][i];
		}
		funktionsværdi[i]=f(A_matrix, b_vektor, Temp, n, m);
	}



	ofstream fil;
	fil.open(fn);

	fil<<"ListPlot[{{"<<t[0]<<","<<funktionsværdi[0]<<"}";
	for(int i=1;i<k+1;i++) fil<<",{"<<t[i]<<","<<funktionsværdi[i]<<"}";
	fil<<"}]";
	fil.close();


}
double f(double A_matrix[nmax][nmax], double b_vektor[nmax], double x[nmax], int n, int m);

double T_InkredsningsAlgoritme(double A_matrix[nmax][nmax], double b_vektor[nmax], double x[nmax], int n, int m, double D, double Epsilon_GS, double s[nmax])
{
	double x1[nmax], x2[nmax], x3[nmax], x4[nmax], a, b; a=0; b=0;

	double LgdDiff, t, k; t=0; k=0;
	double c=(sqrt(5)-1)/2;


	double t1, t2, t3, t4, f1,f2, f3, f4; f1=0; f2=0; f3=0; f4=0;
	t1= 0;t2=0;t3=0;t4=0;
	t2=D;

	for(int i=0;i<m;i++)x1[i] = x[i]+t1*s[i];
	for(int i=0;i<m;i++)x2[i] = x[i]+t2*s[i];
	f1=f(A_matrix,b_vektor,x1,n,m);
	f2=f(A_matrix,b_vektor,x2,n,m);

	while(f2>=f1)
	{
		D=0.1*D;
		t2=t1+D;
		for(int i=0;i<m;i++)x2[i] = x[i]+t2*s[i];
		f2=f(A_matrix,b_vektor,x2,n,m);

	}
	D=2*D;
	t3=t2+D;


	for(int i=0;i<m;i++)x1[i] = x[i]+t1*s[i];
	for(int i=0;i<m;i++)x2[i] = x[i]+t2*s[i];
	for(int i=0;i<m;i++)x3[i] = x[i]+t3*s[i];
	f1=f(A_matrix,b_vektor,x1,n,m);
	f2=f(A_matrix,b_vektor,x2,n,m);
	f3=f(A_matrix,b_vektor,x3,n,m);

	while(f2>f3)
	{
		t1=t2; t2=t3; D=2*D; t3=t2+D;
		for(int i=0;i<m;i++)x2[i] = x[i]+t2*s[i];
		for(int i=0;i<m;i++)x3[i] = x[i]+t3*s[i];
		f2=f(A_matrix,b_vektor,x2,n,m);
		f3=f(A_matrix,b_vektor,x3,n,m);

	}
	a=t1; b=t3;
	t1=0; t2=0;t3=0;t4=0;
	t1=a;
	t4=b;
	t2 = t1+(1-c)*(t4-t1);
	t3 = t4-(1-c)*(t4-t1);

	for(int i=0;i<m;i++)x1[i]=x[i]+t1*s[i];for(int i=0;i<m;i++)x2[i]=x[i]+t2*s[i];
	for(int i=0;i<m;i++)x3[i]=x[i]+t3*s[i];for(int i=0;i<m;i++)x4[i]=x[i]+t4*s[i];

	LgdDiff = fabs(t4-t1);

	f1=0; f2=0; f3=0; f4=0;
	f2=f(A_matrix,b_vektor,x2,n,m);
	f3=f(A_matrix,b_vektor,x3,n,m);


	while(fabs(LgdDiff)>Epsilon_GS)
		{
		if(f2>=f3)
		{
			t1=t2;t2=t3;
			t3=t4-(1-c)*(t4-t1);

			for(int i=0;i<m;i++)x1[i]=x[i]+t1*s[i];for(int i=0;i<m;i++)x2[i]=x[i]+t2*s[i];
			for(int i=0;i<m;i++)x3[i]=x[i]+t3*s[i];for(int i=0;i<m;i++)x4[i]=x[i]+t4*s[i];
			f2=f(A_matrix,b_vektor,x2,n,m);
			f3=f(A_matrix,b_vektor,x3,n,m);

		}
		else
		{
			t4=t3;t3=t2;
			t2=t1+(1-c)*(t4-t1);

			for(int i=0;i<m;i++)x1[i]=x[i]+t1*s[i];for(int i=0;i<m;i++)x2[i]=x[i]+t2*s[i];
			for(int i=0;i<m;i++)x3[i]=x[i]+t3*s[i];for(int i=0;i<m;i++)x4[i]=x[i]+t4*s[i];
			f2=f(A_matrix,b_vektor,x2,n,m);
			f3=f(A_matrix,b_vektor,x3,n,m);

		}
		LgdDiff = fabs(t4-t1);
		}
	t = 0.5*(t4+t1);
	return t;
}
//Kan være fejl med valg af m og/eller n, skiftet rundt eller forkert
double f(double A_matrix[nmax][nmax], double b_vektor[nmax], double x[nmax], int n, int m)
{
	double Ax[nmax], b_Ax[nmax], funktionsværdi;funktionsværdi=0;
	Matrixvektorprodukt_med_m(A_matrix, x, Ax, n,m);
	for(int i=0; i<m;i++) b_Ax[i]=b_vektor[i]-Ax[i];
	for(int i=0;i<m;i++) funktionsværdi += b_Ax[i]*b_Ax[i];
	return funktionsværdi;
}
////2) x* findes ved modificeret Gram-Schmidt
void MinKvad_ved_Mod_GS(double A_matrix[nmax][nmax], double b_vektor[nmax], int n, int m)
{
	double MinKvad[nmax], ATATA[nmax][nmax]; int iii;iii=0;
	double Any[nmax][nmax], Agl[nmax][nmax], r[nmax][nmax], Q[nmax][nmax], QT[nmax][nmax], QTb_vektor[nmax], q[nmax];
	for(int i=0;i< n;i++){for(int j=0;j<m;j++){Any[i][j]=A_matrix[i][j];}}
	for(int i=0;i< n;i++){for(int j=0;j<m;j++){Agl[i][j]=A_matrix[i][j];}}

	for(int k=0;k < m; k++)
	{

		r[k][k]=Søjlelængde_i_matrix(Agl, n, k);
		for(int i=0;i<n;i++) q[i]= Agl[i][k]/r[k][k];

		for(int j=k+1;j<m;j++)r[k][j]=Prikprodukt(q, Agl, n, j);

		for(int j=k+1;j<m;j++)
		{
			for(int i=0;i<n;i++)
			{
				Any[i][j]=Agl[i][j]-r[k][j]*q[i];
			}
		}

		for(int i=0;i< n;i++){for(int j=0;j<m;j++){Agl[i][j]=Any[i][j];}}
		Kopivektorimatrix(q,n, k, Q);
	}
	cout <<endl<<endl<<"Her ses Q matricen: "<<endl;udskriv_matrix(Q, m, m);
	MatrixMatrixprodukt(Q, r, ATATA, n,m);
	cout<<endl<<"Her ses  QR matricen, altså A matricen: "<<endl; udskriv_matrix(ATATA, n, m);
	Transponer(Q, QT, n, m);cout<<endl<<endl<<"Her ses QT matricen: "<<endl;udskriv_matrix(QT,m,n);
	Matrixvektorprodukt_med_m(QT, b_vektor, QTb_vektor, m, n);
	cout<<endl<<endl<<"Her ses b-vektoren:"<<endl; udskriv_vektor(b_vektor, n);
	cout<<endl<<endl<<"her ses QTb vektoren: "<<endl;udskriv_vektor(QTb_vektor, m);
	cout<<endl<<endl<<"Her ses R matricen: "<<endl; udskriv_matrix(r, m ,m);
	Backsub_uden_totalmatrix(r, QTb_vektor, m, MinKvad);
	cout<<"x* ser ud således:"<<endl;udskriv_vektor(MinKvad, m);
}

double Prikprodukt(double qNy[nmax], double Agl[nmax][nmax], int n, int j)
{
	double Prikprodukt; Prikprodukt=0;
	for(int i=0;i<n;i++) Prikprodukt +=qNy[i]*Agl[i][j];
	return Prikprodukt;
}
void Kopivektorimatrix(double vektor[nmax],int n, int j, double Resultatmatrix[nmax][nmax])
{
	for(int i = 0; i<n;i++)
	{
		Resultatmatrix[i][j] = vektor[i];
	}
}
double Søjlelængde_i_matrix(double Matrix[nmax][nmax], int n, int k)
{
	double Længde, Længde2; Længde=0; Længde2=0;
	for(int i=0;i<n;i++) Længde += pow(Matrix[i][k],2);
	Længde2 = sqrt(Længde);
	return Længde2;
	Længde = 0; Længde2=0;
}

double Vektorlængde(double vektor[nmax], int n)
{
	double Længde;Længde=0;
	for(int i=0;i<n;i++) Længde+=pow(vektor[i],2);
	Længde=sqrt(Længde);
	return Længde;
}
////Række-reduktion
void LineærLøsning_uden_totalmatrix(double A_matrix[nmax][nmax], double b_vektor[nmax], int n, double Løsning[nmax])
{
	double Totalmatrix[nmax][nmax];
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n+1; j++)
		{
			if(j < n)
			{
				Totalmatrix[i][j] = A_matrix[i][j];
			}
			else{Totalmatrix[i][n] = b_vektor[i];}
		}
	}
	LineærLøsning(Totalmatrix, Løsning, n);

}
void LineærLøsning(double Totalma[nmax][nmax], double Løsning[nmax], int n)
{
	Gauss(Totalma, n);
	Backsub(Totalma, Løsning, n);

}
void Backsub(double Totalma[nmax][nmax], double x[nmax], int n)
{
	int i, j;
	double sum;
	x[n-1] = Totalma[n-1][n]/Totalma[n-1][n-1];
	for(i = n-2; i >= 0; i--)
	{
		sum = 0;
		for(j = i + 1; j <= n-1; j++)
		{
			sum = sum+Totalma[i][j]*x[j];
		}
		x[i] = (Totalma[i][n] - sum)/Totalma[i][i];

	}
}
void Backsub_uden_totalmatrix(double A_matrix[nmax][nmax], double b_vektor[nmax], int n, double Løsning[nmax])
{
	double Totalmatrix[nmax][nmax];
	for(int i=0;i<n;i++)
	{
		for(int j=0; j<n+1;j++)
		{
			if(j<n){Totalmatrix[i][j]=A_matrix[i][j];}
			else{Totalmatrix[i][n]=b_vektor[i];}
		}
	}
	udskriv_matrix(Totalmatrix, n, n+1);
	Backsub(Totalmatrix, Løsning, n);
}
void DelvisPivotering(double ma[nmax][nmax], int j,  int n)
{
	double max, temp;
	int i, maxpos, k;
	max = fabs(ma[j][j]);
	maxpos = j;
	for(i = j+1; i < n-1; i++)
	{
		if(fabs(ma[i][j]>max))
		{
			maxpos = i;
			max = fabs(ma[i][j]);
		}
	}
	for(k = j; k <= n; k++)
	{
		temp = ma[j][k];
		ma[j][k] = ma[maxpos][k];
		ma[maxpos][k] = temp;
	}
}
void Gauss(double ma[nmax][nmax], int n)
{

	double factor;
	int i, j, k;
	for(j = 0; j <= n-2; j++)
	{
		DelvisPivotering(ma, j, n);
		if(fabs(ma[j][j]) < eps)
		{
			break;
		}
		for(i = j+1; i <= n -1; i++)
		{
			factor = ma[i][j]/ma[j][j];
			ma[i][j] = 0;

			for(k = j+1; k <= n; k++)
			{
			ma[i][k] = ma[i][k] - factor*ma[j][k];
			}
		}
	}
	if(fabs(ma[n-1][n-1]) < eps)
	{

	}
}
void Brugerorientering()
{
	cout << "Velkommen til en besvarelse af eksamen i Numerisk Analyse og Datalogi maj 2020. " << endl;
	cout << "Dette program behandler to hovedemner: Bestemmelse af mindste kvadraters løsning x*, " << endl;
	cout << "til et overbestemt system af lineære ligninger A*x=b, og bestemmelse af egenløsninger " << endl;
	cout << "(lambda_k, v_k), hvor k = 1, 2, ...,n, til en symmetrisk matrix A. " << endl;
	cout << "Frembringelse af datagrundlaget for de forskellige udregninger/metoder kan ske på forskellige måde som du vælger imellem: " << endl;
	cout << "a) Indlæs A, A's rækkeantal n, A's søjleantal m, og vektoren b fra en prækonstrueret datafil" << endl;
	cout << "b) Du indtaster selv n, m elementerne for A og elementerne for b" << endl;
	cout << "c) Simulér målinger (x_i, y_i, z_i), hvor i = 1, 2,...,n, for en multipel regression:" << endl;
	cout << "z_i = a*x_i + b*y_i + c + eps_i, hvor eps_i er simulerede observationer af N(0, sigma)^2 og deraf dannes A og b" << endl;
	cout << "d) Vælg en mindste kvadraters løsning x*, hvor programmet danner A og b således A*x=b er konsistent" << endl << endl;
	cout << "Dernæst beregnes Transpose(A)*A og Transpose(A)*b til brug i nedenstående metoder 1), 4) og i forlængelse 3)" << endl;
	cout << "Bestemmelse af mindste kvadraters løsning x* kan du vælge at skulle udføres ved 3 forskellige metoder:" << endl << endl;
	cout << "1) x* findes som løsning til normalligningerne: Transpose(A)*A*x=Transpose(A)*b" << endl;
	cout << "2) x* findes ved først at udføre QR-faktorisering af A ved modificeret Gram-Schmidt ortogonalisering og dernæst" << endl;
	cout << "løses R*x = Transpose(Q)*b ved backwards substitution" << endl;
	cout << "3) x* findes ved minimering af f(x) = (Det(b-A*x))^2" << endl << endl;
	cout << "Bestemmels af egenløsningerne kan findes ved én metode:" << endl;
	cout << "4) (lambda_k, v_k) af matricen Transpose(A)*A findes ved en udvidet version af potensmetoden" << endl << endl;;
	cout << "Den første del af programmet giver dig mulighed for at vælge imellem de forskellige" << endl;
	cout << "frembringelser af datagrundlag, og den anden del af programmet lader dig vælge" << endl;
	cout << "imellem de 4 algoritmer indenfor de to førnævnte hovedemner" << endl << endl;
	cout << "BEMÆRK: Du skal vælge en specifik forgrening af programmet for at dét kan udføre 3):" << endl;
	cout << "For at udføre 3) skal du vælge 4) først" << endl<< endl;
	cout << "Du har nu muligheden for at vælge imellem de forskellige valgmuligheder af frembringelse af datagrundlag: " << endl;
	cout << "a) Indlæs A, A's rækkeantal n, A's søjleantal m, og vektoren b fra en prækonstrueret datafil" << endl;
	cout << "b) Du indtaster selv n, m elementerne for A og elementerne for b" << endl;
	cout << "c) Simulér målinger (x_i, y_i, z_i), hvor i = 1, 2,...,n, for en multipel regression:" << endl;
	cout << "z_i = a*x_i + b*y_i + c + eps_i, hvor eps_i er simulerede observationer af N(0, sigma)^2 og deraf dannes A og b" << endl;
	cout << "d) Vælg en mindste kvadraters løsning x*, hvor programmet danner A og b således A*x=b er knsistent" << endl << endl;
}
////Datagrundlag d)
void Datagrundlag_d(double A_matrix[nmax][nmax], double b_vektor[nmax], int &n, int &m)
{
	double MinKvad[nmax],Q[nmax][nmax], Q2[nmax][nmax], Q2T[nmax][nmax], Q2Q2T[nmax][nmax], b_p0[nmax], R[nmax][nmax], S[nmax][nmax], sum,bp[nmax], Epsilon[nmax],t[nmax], A_stjerne_matrix[nmax][nmax]; sum=0;
	string fn; fn = "IndlæsAtild.txt";
	cout << "For at skabe en matrix A fra en mindste kvadraters løsning x*, skal du vælge mellem to matricer at indlæse" << endl;
	cout << "Vælg mellem A4 (4*4) og A8 (8*8). For at vælge A4, så indtast 4 og for at vælge A8 så indtast 8: ";cin>> n; cout << endl;
	cout << "Du skal vælge søjleantal m som skal være mindre end " <<n<<". Bemærk at x* kommer til at have m elementer. Indtast m her:";cin>>m;cout<<endl;
	Hent_A_til_d(A_matrix, n, n, fn);
	cout << endl<<"Her ses A matrix som blev indlæst fra en fil"<<endl<<endl; udskriv_matrix(A_matrix,n,n);

	for(int i=0;i<n;i++)
	{
		sum=0;
		for(int j=0;j<n;j++)
		{
		sum+=pow(A_matrix[j][i],2);
		}
		for(int j=0;j<n;j++)S[j][i]=A_matrix[j][i]/sqrt(sum);
	}
	for(int i=0;i<n;i++){for(int j=0;j<m;j++){Q[i][j]=S[i][j];}}

	cout<<endl<<"Her ses Q matricen: "<<endl<<endl; udskriv_matrix(Q,n,m);
	cout <<"Du skal nu indtaste værdier for den øvre trekantsmatrice R." << endl;
	for(int i=0;i<m;i++){for(int j=0;j<m;j++){ if(i<=j){cout<<"Indtast venligst elementet i (række, søjle)= ("<<i+1<<","<<j+1<<"): ";cin>>R[i][j];}}}

	MatrixMatrixprodukt(Q, R, A_stjerne_matrix, n, m);
	cout<< "Her ses A matricen som er udregnet ved A* = QR: "<<endl<<endl;
	udskriv_matrix(A_stjerne_matrix, n,m);

	cout<<endl<<endl<<"Du skal nu indtaste elementer i mindste kvadraters løsningen x*."<<endl;
	for(int i=0;i<m;i++){cout<<"Indtast venligst det "<<i+1<<". element i x*: ";cin>>MinKvad[i];}
	Matrixvektorprodukt(A_stjerne_matrix,MinKvad,bp,n);
	cout<<endl<<endl<<"Her ses bp vektoren:"<<endl;
	udskriv_vektor(bp,n);

	cout<<"Du skal nu indtast elementer for t-vektoren som bruges til at udregne e-vektoren:"<<endl;
	for(int i=0;i <n-m;i++){cout<<"Indtast venligst det "<<i+1<<". element i t: ";cin>>t[i];}
	sum = 0;
	for(int i = 0; i < n;i++)
	{
		for(int j=0; j < n-m; j++)
		{
			Epsilon[i]+= t[j]*A_matrix[i][j+m];
		}
	}
	cout<<endl<<endl<<"Her ses e-vektoren: "<<endl; udskriv_vektor(Epsilon,n);

	for(int i=0;i < n; i++)b_vektor[i]=bp[i]+Epsilon[i];
	cout<<endl<<endl<<"Her ses bp+e = b_vektoren: "<<endl; udskriv_vektor(b_vektor, n);

	cout << "Der foretages en rekvireret designkontrol af c ved at vise projektionen af b på Null(Transpose(A)). Det sker ved formlen"<<endl;
	cout <<"b_p0 = Q2*Transpose(Q2)*b, hvor Q2 er en matrice bestående af de n-m sidste søjler af Q matricen. Først findes Q2: "<< endl;

	for(int i=0; i<n;i++){for(int j=0; j < n-m;j++){Q2[i][j]=S[i][j+m];}}
	cout << "Q2 matricen ses forneden: "<<endl; udskriv_matrix(Q2,n,n-m);
	cout <<"Nu finder vi b's projektion på Null(A): b_p0 = e = Q*Transpose(Q)*b: "<<endl;
	Transponer(Q2, Q2T, n, m);
	MatrixMatrixprodukt(Q2, Q2T,Q2Q2T,n, n);
	cout<<"Matricen Q2*Transpose(Q2) ser ud således: "<<endl; udskriv_matrix(Q2Q2T, n, m);
	Matrixvektorprodukt(Q2Q2T, b_vektor, b_p0, n);
	cout<<"b's profjektion på nulrummet for Transpose(A) ser ud således: "<<endl; udskriv_vektor(b_p0,n);
	cout << endl << "Og fra før havde vi e-vektoren: "<<endl; udskriv_vektor(Epsilon, n);
	cout <<"Som man kan se er det den samme vektor."<<endl;
	for(int i=0; i <n;i++){for(int j=0; j<m;j++){A_matrix[i][j]=A_stjerne_matrix[i][j];}}
}
////Andre hjælpefunktioner
void Valgb(double A[nmax][nmax], double b[nmax], int &n, int &m)
{

	cout << "Indtast venligst antal rækker n i matricen A" << endl; cin >> n;
	cout << "Indtast venligst antal søjler m i matricen A" << endl; cin >> m;

	cout << "Du skal nu indtaste alle elementer i matricen A" << endl << endl;
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			cout << "Du bedes for A indtast  elementet i række nr. " <<i+1<< "og søjle nr. " << j+1<< endl;
			cin >> A[i][j];
		}
	}

	cout << "Du har nu indtastet alle elementer i A. Matricen ser ud således: " << endl << endl;
	udskriv_matrix(A, n, m);


	cout <<endl << "Du skal nu indtaste elementer i vektoren b." << endl;
	for(int i = 0; i < n; i++)
	{
		cout << "Du bedes indtaste det " << i+1 << ". element i vektoren b: " << endl; cin >> b[i];
	}

	cout << "Du har nu indtastet alle elementer i b. vektoren ser ud således: " << endl << endl;
	udskriv_vektor(b, n);

	cout << "Du er nu færdig med at indtaste alle værdier for alle elementer i matricen A og vektoren b" << endl;
}
void randomvektor(double v[nmax], int n)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dis(0.0, 1.0);
	for(int i = 0; i < n; i++) v[i]=dis(gen);
}
void Matrixvektorprodukt(double Matrix[nmax][nmax], double vektor[nmax], double Resultatvektor[nmax], int n)
{
	for(int i = 0; i < n; i++)
		{
			Resultatvektor[i] = 0;
			for(int j = 0; j < n; j++)
			{
				Resultatvektor[i] = Resultatvektor[i] + Matrix[i][j]*vektor[j];
			}
		}
}
void Matrixvektorprodukt_med_m(double Matrix[nmax][nmax], double vektor[nmax], double Resultatvektor[nmax], int n, int m)
{
	for(int i = 0; i < n; i++)
			{
				Resultatvektor[i] = 0;
				for(int j = 0; j < m; j++)
				{
					Resultatvektor[i] = Resultatvektor[i] + Matrix[i][j]*vektor[j];
				}
			}
}
void MatrixMatrixprodukt(double Matrix1[nmax][nmax], double Matrix2[nmax][nmax], double Resultatmatrix[nmax][nmax], int n, int m)
{
	double sum;
		for(int i = 0; i<n; i++)
		{
			for(int j = 0; j<m; j++)
			{
				sum=0;
				for(int k=0; k< n; k++)
				{
					sum+=Matrix1[i][k]*Matrix2[k][j];
					Resultatmatrix[i][j]=sum;
				}
			}
		}

}
void Transponer(double Matrix[nmax][nmax],double Transponeret[nmax][nmax], int n, int m)
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			Transponeret[i][j] = Matrix[j][i];
		}
	}
}
////Metode 4): Bestemmelse af egenløsning
void Egenløsninger(double ATA_matrix[nmax][nmax], int m, double Lambda[nmax], double Egenvektorer[nmax][nmax])
{
	double z[nmax],Epsilon, yGl[nmax], yNy[nmax], Lgl, Lny, LgdyNyMinyGl, sum; Lgl = 0; Lny = 0; sum = 0;
	double Pgl[nmax][nmax], Pny[nmax][nmax];
	int k, N, i, kk; k = 0; i = 0; kk = 0;

	cout << "Du skal nu indtaste tolerancen Epsilon for at finde egenløsning, det skal helst være et meget lille tal. Indtast her: "; cin>> Epsilon;
	cout << endl << "Nu mangler der kun at indtastes det maksimale antal iterationer N, det skal helst være et tilpas stort tal. Indtast her: "; cin>> N;


	for(int i = 0; i < m; i++) for(int j = 0; j < m; j++) Pny[i][j] = ATA_matrix[i][j];
	for(int i = 0; i < m; i++) for(int j = 0; j < m; j++) Pgl[i][j] = Pny[i][j];

	for(k = 0; k < m; k++)
	{
	for(int i = 0; i < m; i++) for(int j = 0; j < m; j++) Pgl[i][j] = Pny[i][j];
	Lgl = 0; Lny = 0; sum = 0;
	randomvektor(z, m);
	NormerVektor(z, yNy, Lny, m);
	do
	{
		Lgl = Lny;
		Kopicib(yNy, m, yGl);
		Matrixvektorprodukt(Pny, yGl, z, m);
		NormerVektor(z, yNy, Lny, m);
		LgdyNyMinyGl = LgdDifVek(yNy, yGl, m);
		kk++;
	}while (fabs(Lny - Lgl) > Epsilon && fabs(LgdyNyMinyGl)> Epsilon	&&	kk < N);
	kk = 0;
	Lambda[k] = Lny;
	for(int ii = 0; ii < m; ii++)sum += pow(yNy[ii], 2);
	for(int ii = 0; ii < m; ii++) Egenvektorer[ii][k] = yNy[ii]/sqrt(sum);
	Pk_kalkulation(Pgl, Egenvektorer, Lny, k, m, Pny);
	}
}
void Pk_kalkulation(double Pgl[nmax][nmax], double Egenvektorer[nmax][nmax], double Lny, int k, int m, double Pny[nmax][nmax])
{
	double v[nmax], v_vT[nmax][nmax];
	for(int i=0; i<m; i++) v[i]=Egenvektorer[i][k];
	for(int i = 0; i < m; i++)for(int j = 0; j < m; j++) v_vT[i][j]=v[i]*v[j];
	for(int i=0;i<m;i++)for(int j=0;j<m;j++)Pny[i][j]= Pgl[i][j] - Lny*v_vT[i][j];
}
double LgdDifVek(double v1[nmax], double v2[nmax], int n)
{
	double DifVek[nmax], sum;
	sum = 0;

	for(int i = 0; i < n; i++)
	{
		DifVek[i] = v1[i] - v2[i];
		sum += pow(DifVek[i], 2);
	}
	sum = pow(sum, 1/2);

	return sum;
}
void NormerVektor(double v[nmax], double vnorm[nmax], double &L, int n)
{
	double maxabs;
	int maxpos;
	maxabs = fabs(v[0]);
	maxpos = 0;

	for(int k = 1; k < n; k++)
	{
		if(fabs(v[k])>maxabs)
		{
			maxabs = fabs(v[k]);
			maxpos = k;
		}
	}
	L = v[maxpos];
	for(int k = 0; k < n; k++)
	{
		vnorm[k] = v[k]/L;
	}
}
void Kopicib(double c[nmax], int n, double b[nmax])
{

	for(int i= 0; i < n; i++)
	{
		b[i] = c[i];
	}
}
//// Datagrundlag c
void Datagrundlag_c_1(double x[nmax],double &deltax,int &mx,double emx[nmax],double y[nmax],double &deltay,int &ny,double eny[nmax],double r[nmax], int &n, int &m)
{
	m = 3;
	cout << "Indtast venligst værdien for x1: "; cin >> x[0];
	cout <<endl<< "Indtast venligst værdien for deltax: "; cin >> deltax;
	cout <<endl<< "Indtast venligst værdien for mx (antal x-værdier): "; cin >> mx;
	cout <<endl<< "Indtast venligst værdien for y1: "; cin >> y[0];
	cout <<endl<< "Indtast venligst værdien for delta y; "; cin >> deltay;
	cout <<endl<< "Indtast venligst værdien for ny (antal y-værdier)"; cin>>ny;
	n = mx*ny;
	Simuler_r(r, n);
	for(int i = 1; i < mx; i++)
	{
		x[i] = x[i-1]+deltax;
		emx[i-1] = 1;
		emx[i] = 1;
	}

	for(int i = 1; i < ny; i++)
	{
		y[i] = y[i-1]+deltay;
		eny[i-1] = 1;
		eny[i] = 1;
	}
}
void Dan_A_og_b(double x[nmax],double deltax,int mx,double emx[nmax],double y[nmax],double deltay,int ny,double eny[nmax],double r[nmax], int n, int m, double A_matrix[nmax][nmax], double b_vektor[nmax])
{
	double Xmatrix[nmax][nmax], Ymatrix[nmax][nmax], sigma, vx[nmax], vy[nmax], a, b, c, Epsilon[nmax];

		cout << "Indtast venligst standardafvigelsen sigma, som skal brugeres til at genere ikke-kontrollérbare tilfældige normalfordelte" << endl;
		cout << "afvigelser ved formlen: sigma*r_i. hvor r_i er de tilfældige obersvationer som blev generet tidligere. Sigma må ikke være for stort" <<endl;
		cout << "Ellers er der for meget støj i den multiple regression. Indtast venligst sigma her: "; cin >> sigma;
		cout <<endl<< "Du skal også indtaste værdierne a, b og c, som skal bruges til den multiple regression for at danne A og b." << endl;
		cout << "Indtast værdien for a: "; cin >>a; cout<<endl<<"Indtast værdien for b: "; cin>>b;cout<<endl<<"Indtast værdien for c";cin>>c;
		cout << endl;
		Ydreprodukt_af_vektorer(eny,x,Xmatrix, ny, mx);
		Ydreprodukt_af_vektorer(y, emx, Ymatrix, ny, mx);
		StakSøjler(vx, Xmatrix, ny, mx);
		StakSøjler(vy, Ymatrix, ny, mx);

		cout<<endl<<"Her ses gitteret af målepunkter for x og y vektoren: "<<endl;

		for(int i=ny-1;i>-1;i--)
		{
			for(int j=-1;j<mx;j++)
			{
				if (j==-1)	cout << y[i]<<"|"<<"\t";
				else
				{
					cout<< "("<<x[j]<<", "<<y[i]<<")\t";
				}
			}
			cout<<endl;
		}

		for(int j=0;j<mx;j++) cout<< "__________"; cout<< endl<<"\t";
		for(int j=0;j<mx;j++) cout<< x[j]<<"   \t";cout<<endl;


		for(int i = 0; i < n; i++)
		{
			A_matrix[i][0] = vx[i];
			A_matrix[i][1] = vy[i];
			A_matrix[i][2] = c;
		}
		for(int i = 0; i < n; i++)
		{
			Epsilon[i] = r[i]*sigma;
			b_vektor[i] = a*vx[i]+b*vy[i]+c+Epsilon[i];
		}
		cout<< endl<<"Her ses r vektoren: "<<endl; udskriv_vektor(r, n);
		cout << endl << endl << "Her ses b_vektoren: " << endl; udskriv_vektor(b_vektor, n);
}
void Simuler_r(double r[nmax], int n)
{
	double xa, xb;
	int i, j;i=0;;xa=0;xb=0;
	j=n;
	do
		{
		xa = (rand()%10001)/10000.0;
		xb =(rand()%10001)/10000.0;
		r[i] = sqrt((-2)*log(xa))*cos(2*M_PI*xb);
		r[j] = sqrt((-2)*log(xa))*sin(2*M_PI*xb);
		i++; j--;

		}while(i < n/2);
}
void NewtonRaphsonFordelingskontrol(double xk[7])
{
    double c, xny, xgl, epsilon = 10e-8, x0[7] = {-1.2, -0.7, -0.4, 0, 0.4, 0.7, 1.2};
    int j = 0, k = 0, i = 0, N =1000;
    cout << "\nNewton-Raphson udført på Phi(xk) = 0.125*k for iterationer af k fra 1 til 7:" << endl;
    do
    {
        k++;
        c = 0.125*k;

        xny = x0[i];

        do
        {
            j++;

            xgl = xny;

            xny = xgl-(FixByTaylor(xgl)-c)/fmrk(xgl);

        } while ((fabs(FixByTaylor(xny)-c) > epsilon) && (j < N));
        cout << "Algortimen konvergerede i nulpunktet x" << k << "= " << xny << endl;
        cout << "Phi(" << xny << ") = " << FixByTaylor(xny) << endl;
        xk[i] = xny;
        i++;
    	}while (k < 7);
}
double FixByTaylor(double x)
{
    double sum = x, term = x, fac = 0, f;
    int k = 1;

    do
    {
        f = sum;

        k = k + 2;

        fac = fac + 1;

        term = term * (-0.5) * ((x*x)/k) * (k-2) * (1/fac);

        sum = sum + term;
    }
    while(f != sum);

    sum = 1/sqrt(2*M_PI)*sum+0.5;

    return sum;
}
void Inddelingafr(double r[nmax], double xk[7], int n)
{
    int Interval[8], Interval_Total; Interval_Total=0;

   for(int i=0;i<8;i++) Interval[i]=0;



    for(int i=0;i<n;i++)
    {
    	if(r[i]<xk[0]) Interval[0]++;

    }


for(int k=1; k<6;k++)
{
	for(int i=0;i<n;i++)
	{
		if(r[i]> xk[k-1] && r[i]<=xk[k]) Interval[k]++;
	}
}

for(int i=0;i<n;i++)
{
	if(r[i]>xk[6]) Interval[7]++;
}

    cout << "\nVærdierne af de forskellige r er inddelt i intervallerne således: " << endl;

for(int i=0;i< 8;i++)
{
	cout<<"Antal i interval nr. "<<i+1<< " er: "<<Interval[i]<<endl;
	Interval_Total+=Interval[i];

}


    cout << "Total: " << Interval_Total << endl;
}
void Ydreprodukt_af_vektorer(double vektor_a[nmax], double vektor_b[nmax], double Resultatmatrix[nmax][nmax], int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			Resultatmatrix[i][j] = vektor_a[i]*vektor_b[j];
		}
	}
}
void StakSøjler(double vektor[nmax], double Matrix[nmax][nmax], int n, int m)
{
	int a;

	for(int j = 0; j < m; j++)
	{
		a = j*n;
		for(int i = 0; i < n; i++)
		{
			vektor[a+i] = Matrix[i][j];
		}
	}
}
double fmrk(double x)
{
	double resultat;
	resultat = (1/(sqrt(2*M_PI)))*exp(-1/2*pow(x,2)); // Her defineres differentialkvotienten til f
	return resultat;
}
double Simpson(double a, double b, int ndelInt)
{
	double h, sum1, sum2 = 0, sum3 = 0, I = 0;

	h = (b-a)/ndelInt;

	sum1 = fmrk(a) + fmrk(b);

	for(int i = 1; i <= ndelInt-1; i++)
	{
		sum2 = sum2+fmrk(a+i*h);
	}
	for(int k=1; k <= ndelInt; k++)
	{
		sum3 = sum3 + fmrk(a - (h/2) + k * h);
	}
	I = (h/6)*(sum1 + 2*sum2 + 4*sum3); // I = Bestemt integral tilnærmelse ved Simpson
	return I;
}
////Udskrivningsfunktioner
void udskriv_matrix(double A[nmax][nmax], int n, int m) {
	for(int i =0; i < n; i++) {
		for (int j = 0; j< m; j++)
		{
			if(fabs(A[i][j]) > 0.0001)
			{
			cout <<  setw(15) << A[i][j];
			}
			else
			{
				cout <<  setw(15) << "0";
			}
		}
		cout << endl;
	}
}
void udskriv_Totalmatrix(double A[nmax][nmax+1], int n, int m)
{
	for(int i =0; i < n; i++) {
		for (int j = 0; j< m+1; j++) {
			cout <<  setw(15) << A[i][j];
		}
		cout << endl;
	}
}
void udskriv_vektor(double v[nmax], int n)
{
	for(int i = 0; i <n; i++)
	{
		cout<< v[i] << endl;
	}
}
////fstream funktioner
void hent_Aogb(double A[nmax][nmax], double b[nmax],int &n, int &m, string fn)
{
	ifstream fil;
	fil.open(fn);
	fil>>n;
	fil>>m;
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<m+1;j++)
		{
			if(j<m)
			{
			fil>>A[i][j];
			}
			else
			{
			fil>>b[i];
			}
		}
	}
	fil.close();
}
void Hent_A_til_d(double A[nmax][nmax], int n, int m, string fn)
{
	double temp[nmax][nmax];
	ifstream fil;
	fil.open(fn);
		for (int i=0;i<n;i++)
		{
			if(n == 4)
			{
			for(int j=0;j<m;j++)fil>>A[i][j];
			for(int j = 0; j < m; j++) fil>>temp[i][j];
			}
			else
			{
				for(int j=0; j< m;j++)fil>>A[i][j];
			}
		}
		fil.close();
}
void gem_Aogb(double A[nmax][nmax], double b_vektor[nmax],int n, int m, string fn)


{
	ofstream fil;
	fil.open(fn);
	fil<<n<<endl;
	fil<<m<<endl;
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<m+1;j++)
		{
			if(j<m)
			{
				fil<<A[i][j]<<"\t";
			}
			else
			{
				fil<<b_vektor[i];
			}
		}
		fil<<endl;
	}
	fil.close();
}
void hent_matrix(double M[nmax][nmax],int &n, string fn){
	ifstream fil;
	fil.open(fn);
	fil>>n;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++) fil>>M[i][j];
	}
	fil.close();
}
int vaelg_A(double A[nmax][nmax]) {
	string afil;
	int n;
	do {
		cout << "Indtast n: ";
		cin >> n;
		switch (n) {
		case 4:
			afil = "a4.txt";
			break;
		case 8:
			afil = "a8.txt";
			break;
		case 16:
			afil = "a16.txt";
			break;
		default:
			cout << "Indtast 4, 8 eller 16.\n";
			continue;
		}
	} while (afil == "");
	hent_matrix(A, n, afil);
	return n;
}
