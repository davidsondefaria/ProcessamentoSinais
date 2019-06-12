#include <stdio.h>
#include <math.h>

#define pi 4*atan(1)
#define TAM 257		//quantidade de pontos do rex
#define NF 500		//numero do filtro: 20, 100, 500

main(){
	int i, j, k, N;
	N = 4*(TAM-1);
	int nm1,nd2,m,l,le,inp,le2,jm1, ind;
	double tr,ti,ur,si,sr,ui;	
	double fx[TAM], fy[TAM];
	double rex[N], imx[N], temp[N];
	
	//arquivo de leitura do filtro
	FILE *farq;		
	farq = fopen("pontoAjustado.dat", "r");
	for(i=0; i<TAM; i++){
		fscanf(farq, "%lf\t%lf\n", &fx[i], &fy[i]);
	}
	fclose(farq);
	
	//lê rex
	for(i=0; i<N; i++){
		rex[i] = 0;
		imx[i] = 0;
	}
	for(i=0; i<TAM; i++){
		rex[i] = fy[i];
	}

//FFT INVERSA
	//inverte sinal da imx
	for(k=0; k<N; k++)
		imx[k] = -imx[k];
	
	//calcula FFT
	nm1 = N-1;
	nd2 = N/2;
	m = log(N)/log(2);
	j = nd2;
	k=0;
	
	//bit reverse
	for(i=1; i<N-1; i++){
		if(i < j){
			tr = rex[j];
			ti = imx[j];
			rex[j] = rex[i];
			imx[j] = imx[i];
			rex[i] = tr;
			imx[i] = ti;
		}
		k = nd2;
		while(k<=j){
			j = j-k;
			k = k/2;
		}
		j = j+k;
	}
	//loop dos estagios
	for(l=1; l<=m; l++){
		le = pow(2, l);
		le2 = le/2;
		ur = 1.0;
		ui = 0.0;
		sr = cos(pi/le2);		//calcula valores de senos e cossenos
		si = -sin(pi/le2);
		for(j=1; j<=le2; j++){	//calcula cada sub dft
			jm1 = j-1;
			for(i=jm1; i<=nm1; i+=le){				//loop para cada butterfly
				inp = i+le2;
				tr = rex[inp]*ur - imx[inp]*ui;		//calculo da butterfly
				ti = rex[inp]*ui + imx[inp]*ur;
				rex[inp] = rex[i] - tr;
				imx[inp] = imx[i] - ti;
				rex[i] = rex[i] + tr;
				imx[i] = imx[i] + ti;

			}
			tr = ur;
			ur = tr*sr - ui*si;
			ui = tr*si + ui*sr;
		}
	}

	//divide dominio por N e muda sinal de imx
	for (i=0;i<N;i++){
		rex[i] = rex[i]/(double)N;
		imx[i] = -imx[i]/(double)N;
	}
	
	//impressao teste Impulso
	FILE *imparq;		
	imparq = fopen("Impulso500.dat", "w");
	for(i=0; i<N/2; i++){
		fprintf(imparq, "%d\t%f\n", i,  sqrt(rex[i]*rex[i] + imx[i]*imx[i]));
	}
	fclose(imparq);	
	
//Janelamento
	for(i=0; i<N; i++){
		ind = i+NF/2;
		if(ind > N-1)
			ind = ind - N;
		temp[ind] = rex[i];
		rex[i] = 0;
		imx[i] = 0;
	}
	for(i=0; i<N; i++)
		rex[i] = temp[i];

	for(i=0; i<N; i++){
		if (i<=NF){
			rex[i] = rex[i]*(0.54 - 0.46*cos((2*pi*i)/NF));
		}
		else{
			rex[i] = 0;
			imx[i] = 0;
		}
	}
	FILE *rexarq;		
	rexarq = fopen("rexJanelamento500.dat", "w");
	for(i=0; i<N/2; i++){
		fprintf(rexarq, "%d\t%f\n", i, sqrt(rex[i]*rex[i]+imx[i]*imx[i]));
	}
	fclose(rexarq);	

//calcula FFT
	nm1 = N-1;
	nd2 = N/2;
	m = log(N)/log(2);
	//printf("N: %d, m: %d", N, m);
	j = nd2;
	k=0;
	//bit reverse
	for(i=1; i<N-1; i++){
		if(i < j){
			tr = rex[j];
			ti = imx[j];
			rex[j] = rex[i];
			imx[j] = imx[i];
			rex[i] = tr;
			imx[i] = ti;
		}
		k = nd2;
		while(k<=j){
			j = j-k;
			k = k/2;
		}
		j = j+k;
	}
	//loop dos estagios
	for(l=1; l<=m; l++){
		le = pow(2, l);
		le2 = le/2;
		//printf("%d\n", le2);
		ur = 1.0;
		ui = 0.0;
		sr = cos(pi/le2);		//calcula valores de senos e cossenos
		si = -sin(pi/le2);
		for(j=1; j<=le2; j++){	//calcula cada sub dft
			jm1 = j-1;			
			for(i=jm1; i<=nm1; i+=le){				//loop para cada butterfly
				inp = i+le2;
				tr = rex[inp]*ur - imx[inp]*ui;		//calculo da butterfly
				ti = rex[inp]*ui + imx[inp]*ur;
				rex[inp] = rex[i] - tr;
				imx[inp] = imx[i] - ti;
				rex[i] = rex[i] + tr;
				imx[i] = imx[i] + ti;
			}
			tr = ur;
			ur = tr*sr - ui*si;
			ui = tr*si + ui*sr;
		}
	}
	FILE *fftarq;		
	fftarq = fopen("FFT500.dat", "w");
	for(i=0; i<(N/4)+1; i++){
		fprintf(fftarq, "%f\t%f\n",((float) i/(N/2)), 2*sqrt(rex[i]*rex[i]+imx[i]*imx[i]));
	}
	fclose(fftarq);	
}
















