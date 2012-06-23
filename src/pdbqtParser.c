#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <regex.h>

#define MAX_ATOMS 5000
#define ALFA 2
#define SQR(x) x*x
#define length(x) (sizeof(x)/sizeof(x[0]))

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)       __builtin_expect((x),0)

#define MAX(a, b) (a > b) ? a : b
#define MIN(a, b) (a < b) ? a : b
//#define MAX(a, b) (a ^ ((a ^ b) & -(a < b)))
//#define MIN(a, b) (b ^ ((a ^ b) & -(a < b)))

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
     float x;
     float y;
     float z;
     float q;
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

inline point getCenter (atom* molecule) {
	unsigned int i;
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

int readPdbqt(char* fname, molecule* new_molecule) {

	int i = 0;
	char line[512];

	new_molecule->natoms = 0;

/*
	if( ) ) {
		fprintf(stderr, "Not enougth memory avaliable");
	}
*/
	FILE *file = fopen(fname,"r");
	if (!file) {
		printf("Error opening output file\n");
		return 1;
	}

	while (fgets(line,sizeof(line),file)){
		if (match(line,"ATOM") || match(line,"HETATM")) {
			sscanf(&line[31], "%f", &(new_molecule->atoms[i]).x);
			sscanf(&line[38], "%f", &(new_molecule->atoms[i]).y);
			sscanf(&line[46], "%f", &(new_molecule->atoms[i]).z);
			sscanf(&line[70], "%f", &(new_molecule->atoms[i]).q);
			sscanf(&line[77], "%s", (new_molecule->atoms[i]).type);
			i++;
		}
	}
	fclose(file);
	
	(new_molecule->natoms) = i;

	return 0;
}

#define DISTANCE_SQR(A, B) ((pow((A->atoms[i]).x - (B->atoms[i]).x, 2)) \
                            + (pow((A->atoms[i]).y - (B->atoms[i]).y, 2)) \
                            + (pow((A->atoms[i]).z - (B->atoms[i]).z, 2)))

double computeKernel(molecule* m1, molecule* m2) {

    unsigned int i, j;
	double kernelValue=0;

    double norm1=0; 
    double norm2=0;
	double alfa = 2;

	for (i=0; i<(m1->natoms); i++) {
		for (j=0; j<(m2->natoms); j++) {
         
			kernelValue +=  (m1->atoms[i]).q 
                            * (m2->atoms[j]).q 
                            * exp((-alfa/2) * DISTANCE_SQR(m1, m2));

		}
	}

	for (i=0; i<(m1->natoms); i++) {
		norm1 += SQR((m1->atoms[i]).q);
	}
	norm1 = sqrt(norm1);

	for (j=0; j<(m2->natoms); j++) {
		norm2 += SQR((m2->atoms[j]).q);
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

    if(argc != 3) {
        fprintf(stderr, "./pdqtParser molecule_A.pdbqt molecule_B.pdbqt\n");
        return 1;
    }

    unsigned int i;
    char* fname1 	= argv[1];
    char* fname2 	= argv[2];

	molecule molec1;
    readPdbqt(fname1, &molec1);

	molecule molec2;
    readPdbqt(fname2, &molec2);


	for (i=0; i<molec1.natoms; i++) {
		printf("%d\t%lf\t%lf\t%lf\t%lf\t%s\n",i,
                molec1.atoms[i].x,molec1.atoms[i].y,molec1.atoms[i].z,
                molec1.atoms[i].q,molec1.atoms[i].type);
	}

    for (i=0; i<molec2.natoms; i++) {
		printf("%d\t%lf\t%lf\t%lf\t%lf\t%s\n",i,
                molec2.atoms[i].x,molec2.atoms[i].y,molec2.atoms[i].z,
                molec2.atoms[i].q,molec2.atoms[i].type);
	}

	double kernel = computeKernel(&molec1, &molec2);
	printf("Kernel: %lf\n",kernel);

//	computeGrid(molec1);
//	computeGrid(molec2);


	return 0;
}


