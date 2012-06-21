#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <regex.h>

#define MAX_ATOMS 5000
#define ALFA 2
#define SQR(x) x*x
#define length(x) (sizeof(x)/sizeof(x[0]))

int match(const char *string, char *pattern) {
	int status;
	regex_t re;

	if (regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0) { return(0); } // Report error. 
	status = regexec(&re, string, (size_t) 0, NULL, 0);
	regfree(&re);
	if (status != 0) { return(0); }      // Report error.
	return(1);
}

typedef struct {
     double x;
     double y;
     double z;
     double q;
     char type[2];
     char dummy;
} atom;

typedef struct {
	unsigned int natoms;
	atom atoms[MAX_ATOMS];
} molecule;

typedef struct {
	float x;
	float y;
	float z;
} point;

point getCenter (atom* molecule) {
	int i;
	point center;
	center.x=0;
	center.y=0;
	center.z=0;
	for (i=0; i<sizeof(atom); i++) {
		center.x+=molecule[i].x;
		center.y+=molecule[i].y;
		center.z+=molecule[i].z;
	}
	center.x/=sizeof(molecule);
	center.y/=sizeof(molecule);
	center.z/=sizeof(molecule);
	return center;
}

molecule readPdbqt(char* fname) {

	int i=0;
	char line[512];

	printf("Hi");
	molecule new_molecule;
	new_molecule.natoms = 0;

/*
	if( ) ) {
		fprintf(stderr, "Not enougth memory avaliable");
	}
*/
	FILE *file = fopen(fname,"r");
	if (!file) {
		printf("Error opening output file\n");
		return new_molecule;
	}

	printf("Hi");
	while (fgets(line,sizeof(line),file)){
		if (match(line,"ATOM") || match(line,"HETATM")) {
			sscanf(&line[31],"%lf",new_molecule.atoms[i].x);
			sscanf(&line[38],"%lf",new_molecule.atoms[i].y);
			sscanf(&line[46],"%lf",new_molecule.atoms[i].z);
			sscanf(&line[70],"%lf",new_molecule.atoms[i].q);
			sscanf(&line[77],"%s",new_molecule.atoms[i].type);
			i++;
		}
	}
	fclose(file);
	
	new_molecule.natoms = i;

	return new_molecule;
}

double computeKernel(atom* m1, atom* m2) {
	int i,j;
	double kernelValue=0;
	double norm1=0;
	double norm2=0;
	double alfa = 2.;
	for (i=0; i<sizeof(*m1); i++) {
		for (j=0; j<sizeof(*m2); j++) {
			kernelValue += m1[i].q * m2[j].q * exp((-alfa/2) * ( pow((m1[i].x-m2[j].x),2) + pow(m1[i].y-m2[j].y,2) + pow(m1[i].z-m2[j].z,2)));
		}
	}
	for (i=0; i<sizeof(*m1); i++) {
		norm1 += SQR(m1[i].q);
	}	
	norm1 = sqrt(norm1);

	for (j=0; j<sizeof(*m2); j++) {
		norm2 += SQR(m2[j].q);
	}	
	norm2 = sqrt(norm2);

	printf("%lf\t%lf\t%lf\n",kernelValue,norm1,norm2);
	kernelValue = kernelValue/(norm1*norm2);
	
	return kernelValue;
}

void computeGrid(atom* m) {

	point center = getCenter(m);

	printf("%lf\t%lf\t%lf\n",center.x,center.y, center.z);
	


}

int main (int argc, char* argv[]) {

	char* fname1 	= argv[1];
	char* fname2 	= argv[2];

	molecule molec2 = readPdbqt(fname2);
	molecule molec1 = readPdbqt(fname1);
	printf("bu");
	int i;
	for (i=0; i<molec2.natoms; i++) {

		printf("%d\t%lf\t%lf\t%lf\t%lf\t%s\n",i,molec2.atoms[i].x,molec2.atoms[i].y,molec2.atoms[i].z,molec2.atoms[i].q,molec2.atoms[i].type);

	}
	

//	double kernel = computeKernel(molec1, molec2);
//	printf("Kernel: %lf\n",kernel);

//	computeGrid(molec1);
//	computeGrid(molec2);


	return 0;
}


