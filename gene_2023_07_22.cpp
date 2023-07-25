#include <iostream>
#include <time.h>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

const int gPopulationCount = 80; // �ڸs�ƶq
const int gNumberOfGenes = 8;
const int gExpansionCount = 10000; // �X�R�ƶq�C���ܼ����u�X�{�b��L�ܼƫŧi�̡C
int gExternalCount = 0; // ��e�M���ܮɬ����B�~�W�[���ƶq

double Generate(double min, double max) {
    // �ͦ� min �� max �������H����
    double c;
    double Beta;
    Beta = (double) rand() / (RAND_MAX + 1.0);
    c = min + Beta * (max - min);
    c = round(c * 1000000) / 1000000.0; // �u�O�d�e6��A�᭱�˱�C
    return c;
}

double GetUpperLimit(int gene) {
	// ��]�W��
	switch (gene) {
        case 0:
        case 1:
        case 2:
            return 10000;
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            return 1000;
        default:
            return 0;
    }
}

double GetLowerLimit(int gene) { 
	// ��]�U��
    switch (gene) {
        case 0:
            return 100;
        case 1:
        case 2:  
			return 1000;
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            return 10;
        default:
            return 0;
    }
}

double Judge(double g1, double g2, double g3, double g4, double g5, double g6, double g7, double g8) {
    // �p��P�_�����~�`�Z��
    double sum = 0;
    double temp;
    if ((temp = 1 - 0.0025 * (g4 + g6)) < 0) {
        sum += abs(temp);
    }
    if ((temp = 1 - 0.0025 * (g5 + g7 - g4)) < 0) {
        sum += abs(temp);
    }
    if ((temp = 1 - 0.01 * (g8 - g5)) < 0) {
        sum += abs(temp);
    }
    if ((temp = g1 * g6 - 833.33252 * g4 - 100 * g1 + 83333.333) < 0) {
        sum += abs(temp);
    }
    if ((temp = g2 * g7 - 1250 * g5 - g2 * g4 + 1250 * g4) < 0) {
        sum += abs(temp);
    }
    if ((temp = g3 * g8 - 1250000 - g3 * g5 + 2500 * g5) < 0) {
        sum += abs(temp);
    }
    return sum;
}

void GenerateGene(double population[][gNumberOfGenes], double fitnessPunishment[]) {
    double g1, g2, g3, g4, g5, g6, g7, g8, judgement;
    for (int i = 0; i < gPopulationCount + gExternalCount; i++) {
        g1 = Generate(GetLowerLimit(0), GetUpperLimit(0));
        g2 = Generate(GetLowerLimit(1), GetUpperLimit(1));
        g3 = Generate(GetLowerLimit(2), GetUpperLimit(2));
        g4 = Generate(GetLowerLimit(3), GetUpperLimit(3));
        g5 = Generate(GetLowerLimit(4), GetUpperLimit(4));
        g6 = Generate(GetLowerLimit(5), GetUpperLimit(5));
        g7 = Generate(GetLowerLimit(6), GetUpperLimit(6));
        g8 = Generate(GetLowerLimit(7), GetUpperLimit(7));
		
		if(Judge(g1, g2, g3, g4, g5, g6, g7, g8) != 0) {
			i--;
			continue;
		}

        // �p��P�_�����~�`�Z�� 
        judgement = Judge(g1, g2, g3, g4, g5, g6, g7, g8);
        fitnessPunishment[i] = judgement * 10000000000.0; // �Y���~�t�A�N��~�t�ȩ�j

        // �s�J�}�C 
        population[i][0] = g1; 
        population[i][1] = g2; 
        population[i][2] = g3; 
        population[i][3] = g4; 
        population[i][4] = g5; 
        population[i][5] = g6; 
        population[i][6] = g7; 
        population[i][7] = g8;
    }
	// ��X�ثe�ڸs 
	cout << "��N��]: \n";
	for(int j = 0 ; j < gPopulationCount + gExternalCount ; j++){
		cout << "G " << j << " : ";
		for(int k = 0 ; k < gNumberOfGenes ; k++){
            cout << population[j][k] << "    " ;  
    	}
    	cout<<endl;
    }
}

void Sortpopulation(double population[][gNumberOfGenes], double fitnessValues[], double fitnessPunishment[]) {
	int totalPopulationSize = gPopulationCount + gExternalCount;
	double sortedFitnessValues[gPopulationCount + gExpansionCount];    // �x�s�Ƨǫ᪺�}�C
	int sortedIndices[gPopulationCount + gExpansionCount] = {0};       // �x�s�Ƨǫ�ҹ����}�C��m
	
	for(int i = 0; i < gPopulationCount + gExternalCount ; i++){
		fitnessValues[i] = population[i][0] +  population[i][1] + population[i][2] + fitnessPunishment[i];  
	}
	
	memcpy(sortedFitnessValues, fitnessValues, sizeof(double) * totalPopulationSize);
	sort(sortedFitnessValues, sortedFitnessValues + totalPopulationSize);
	// �O���ҹ������ڸs��m
	for(int i = 0 ; i < totalPopulationSize ; i++){
		for(int j = 0; j < gPopulationCount + gExternalCount ; j++){
			if(sortedFitnessValues[i] == fitnessValues[j]){
				sortedIndices[i] = j;
			}
		}
	} 
	//�N�ڸs�Ƨǫ᪺��m��s
	double populationTemp[gPopulationCount + gExpansionCount][gNumberOfGenes] = {0}; // �Ȧs�Ƨǫ�ڸs�}�C
	for(int i = 0 ; i < totalPopulationSize ; i++){
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			int temp = sortedIndices[i];
			populationTemp[i][j] = population[temp][j];
		}
	}
	//�N��m��s�᪺�Ȧs�}�C ���s��^��l�ڸs�}�C 
	for(int i = 0 ; i < totalPopulationSize ; i++){
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			population[i][j] = populationTemp[i][j];
		}
	} 
	//�N�Ƨǫ᪺�A����ƭȡA���s��^��l�A����ƭȰ}�C 
	memcpy(fitnessValues, sortedFitnessValues, sizeof(double) * totalPopulationSize);
}

double CalculateFitnessTotal(double fitnessValues[]) {
	double sum = 0;
    for (int i = 0; i < gPopulationCount; i++) {
        sum += fitnessValues[i];
    }
	return sum;
}

void CalculateChosenProbability(double chosenProbability[], double fitnessValues[], double fitnessTotalValues) {
    for (int i = 0; i < gPopulationCount; i++) {
		if(fitnessTotalValues == 0)
			chosenProbability[i] = 0;
		else
			chosenProbability[i] = fitnessValues[i] / fitnessTotalValues;
    }
}

void CalculateAccumulatedProbabilities(double accumulatedProbabilities[], double chosenProbability[]) {
    double temp = 0;
    for (int i = 0; i < gPopulationCount; i++) {
        temp += chosenProbability[i];
        accumulatedProbabilities[i] = temp;
    }
}

int SelectIndividualForCrossover(double accumulatedProbabilities[]) {
    int selectedForCrossover = 0;
	double randNum;
    randNum = (double)rand() / (RAND_MAX + 1.0);
    for (int i = 1; i < gPopulationCount ; i++) {
        if (randNum > accumulatedProbabilities[i - 1] && randNum <= accumulatedProbabilities[i]) {
            selectedForCrossover = i;
            selectedForCrossover = gPopulationCount - selectedForCrossover + 1; // ��p����
            break;
        }
    }
    return selectedForCrossover;
}

void AddIndividualToPopulation(double population[][gNumberOfGenes], int gene) {
    gExternalCount++;
    for (int i = 0; i < gNumberOfGenes; i++) {
        population[gPopulationCount + gExternalCount - 1][i] = population[gene][i];
    }
}

double singlePointInnerProductX(double B, double x1, double x2) {
    // ���I���n���
    double ret;
    ret = B * x1 + (1 - B) * x2;
    ret = (int)(ret * 1000000) / 1000000.0;
    return ret;
}

double singlePointInnerProductY(double B, double y1, double y2) {
    // ���I���n��� (Y��)
    double ret;
    ret = (1 - B) * y1 + B * y2;
    ret = (int)(ret * 1000000) / 1000000.0;
    return ret;
}

void Crossover(double population[][gNumberOfGenes], double chosenProbability[], double accumulatedProbabilities[]) {
	int crossRate = static_cast<int>(gPopulationCount * 0.8);
	int selectedForCrossover1, selectedForCrossover2; // �Q��ܰ���e������
	int selectedCrossGene; // �Q��ܰ���e����]
	double gene1, gene2; // ��e����]
	for(int i = 0 ; i < crossRate ; i++) {
		// ����ֿn���v�A��X����e������
		selectedForCrossover1 = SelectIndividualForCrossover(accumulatedProbabilities);
		selectedForCrossover2 = SelectIndividualForCrossover(accumulatedProbabilities);
    	
    	selectedCrossGene = rand() % gNumberOfGenes;
    	gene1 = population[selectedForCrossover1][selectedCrossGene];
    	gene2 = population[selectedForCrossover2][selectedCrossGene];
		
		// �n��e�����鱵�b�ڸs����
        AddIndividualToPopulation(population, selectedForCrossover1);
        AddIndividualToPopulation(population, selectedForCrossover2);
		
		double Beta = static_cast<double>(rand()) / (RAND_MAX + 1.0); // ���� 0 ~ 1 �������B�I��
		// ���n(���I���n)  (100+2-2=100)
		population[gNumberOfGenes + gExternalCount - 2][selectedCrossGene] = singlePointInnerProductX(Beta, gene1, gene2);
		// ���n(���I���n)  (100+2-1=101)
    	population[gNumberOfGenes + gExternalCount - 1][selectedCrossGene] = singlePointInnerProductY(Beta, gene1, gene2);
		// ���I��Ҧ���������
		for (int j = selectedCrossGene; j < gNumberOfGenes; j++) {
			double temp = population[selectedForCrossover1][j];
			population[selectedForCrossover1][j] = population[selectedForCrossover2][j];
			population[selectedForCrossover2][j] = temp;
		}
	}
}

void Tago(double population[][gNumberOfGenes], double accumulatedProbabilities[]) {
	// ���Фf��u�A�Ω��]�ƥص��� 8 ��
	int sizeofTago = 16;
    double tagArr[sizeofTago][gNumberOfGenes] = {0};
	tagArr[0][0]=0;tagArr[0][1]=0;tagArr[0][2]=0;tagArr[0][3]=0;tagArr[0][4]=0;tagArr[0][5]=0;tagArr[0][6]=0;tagArr[0][7]=0;
	tagArr[1][0]=0;tagArr[1][1]=0;tagArr[1][2]=0;tagArr[1][3]=0;tagArr[1][4]=0;tagArr[1][5]=0;tagArr[1][6]=0;tagArr[1][7]=1;
	tagArr[2][0]=0;tagArr[2][1]=0;tagArr[2][2]=0;tagArr[2][3]=1;tagArr[2][4]=1;tagArr[2][5]=1;tagArr[2][6]=1;tagArr[2][7]=0;
	tagArr[3][0]=0;tagArr[3][1]=0;tagArr[3][2]=0;tagArr[3][3]=1;tagArr[3][4]=1;tagArr[3][5]=1;tagArr[3][6]=1;tagArr[3][7]=1;
	tagArr[4][0]=0;tagArr[4][1]=1;tagArr[4][2]=1;tagArr[4][3]=0;tagArr[4][4]=0;tagArr[4][5]=1;tagArr[4][6]=1;tagArr[4][7]=0;
	tagArr[5][0]=0;tagArr[5][1]=1;tagArr[5][2]=1;tagArr[5][3]=0;tagArr[5][4]=0;tagArr[5][5]=1;tagArr[5][6]=1;tagArr[5][7]=1;
	tagArr[6][0]=0;tagArr[6][1]=1;tagArr[6][2]=1;tagArr[6][3]=1;tagArr[6][4]=1;tagArr[6][5]=0;tagArr[6][6]=0;tagArr[6][7]=0;
	tagArr[7][0]=0;tagArr[7][1]=1;tagArr[7][2]=1;tagArr[7][3]=1;tagArr[7][4]=1;tagArr[7][5]=0;tagArr[7][6]=0;tagArr[7][7]=1;
	tagArr[8][0]=1;tagArr[8][1]=0;tagArr[8][2]=1;tagArr[8][3]=0;tagArr[8][4]=1;tagArr[8][5]=0;tagArr[8][6]=1;tagArr[8][7]=0;
	tagArr[9][0]=1;tagArr[9][1]=0;tagArr[9][2]=1;tagArr[9][3]=0;tagArr[9][4]=1;tagArr[9][5]=0;tagArr[9][6]=1;tagArr[9][7]=1;
	tagArr[10][0]=1;tagArr[10][1]=0;tagArr[10][2]=1;tagArr[10][3]=1;tagArr[10][4]=0;tagArr[10][5]=1;tagArr[10][6]=0;tagArr[10][7]=0;
	tagArr[11][0]=1;tagArr[11][1]=0;tagArr[11][2]=1;tagArr[11][3]=1;tagArr[11][4]=0;tagArr[11][5]=1;tagArr[11][6]=0;tagArr[11][7]=1;
	tagArr[12][0]=1;tagArr[12][1]=1;tagArr[12][2]=0;tagArr[12][3]=0;tagArr[12][4]=1;tagArr[12][5]=1;tagArr[12][6]=0;tagArr[12][7]=0;
	tagArr[13][0]=1;tagArr[13][1]=1;tagArr[13][2]=0;tagArr[13][3]=0;tagArr[13][4]=1;tagArr[13][5]=1;tagArr[13][6]=0;tagArr[13][7]=1;
	tagArr[14][0]=1;tagArr[14][1]=1;tagArr[14][2]=0;tagArr[14][3]=1;tagArr[14][4]=0;tagArr[14][5]=0;tagArr[14][6]=1;tagArr[14][7]=0;
	tagArr[15][0]=1;tagArr[15][1]=1;tagArr[15][2]=0;tagArr[15][3]=1;tagArr[15][4]=0;tagArr[15][5]=0;tagArr[15][6]=1;tagArr[15][7]=1;
	
	int crossRate = static_cast<int>(gPopulationCount * 0.8);
	int selectedForCrossover1, selectedForCrossover2; // �Q��ܰ���e������
	double g1, g2, g3, g4, g5, g6, g7, g8, judgement;
	double tagoGeneArr[2][gNumberOfGenes]; // �Фf���L
	double tagoPopulationMatrix[sizeofTago][gNumberOfGenes] = {0}; // �Фf�}�C 
	double tagoFitnessValue[sizeofTago] = {0}; // �Фf�A����ƭ�
	double tagoFitnesscPow[sizeofTago] = {0}; // �Фf FitnessValue ��������@ 
	double ef1[gNumberOfGenes] = {0};
	double ef2[gNumberOfGenes] = {0};
	int optimalLevel[gNumberOfGenes] = {0};
	double tagoFinalGene[gNumberOfGenes] = {0}; // �Фf�̫��]��
	for(int i = 0 ; i < crossRate ; i++) {
		// ����ֿn���v�A��X����e������
		selectedForCrossover1 = SelectIndividualForCrossover(accumulatedProbabilities);
		selectedForCrossover2 = SelectIndividualForCrossover(accumulatedProbabilities);
		// �ҿ� 2 �զs�J�}�C
		for(int j = 0 ; j < gNumberOfGenes ; j++) {
        	tagoGeneArr[0][j] = population[selectedForCrossover1][j];
			tagoGeneArr[1][j] = population[selectedForCrossover2][j];
        }
		//�N��]�ȱa�J�Фf��
    	for(int j = 0 ; j < sizeofTago ; j++) {
    		for(int k = 0 ; k < gNumberOfGenes ; k++){
    			int temp = tagArr[j][k]; 
    			tagoPopulationMatrix[j][k] = tagoGeneArr[temp][k];
			}
		}
		// �p�� fitness value (�A����ƭ�)
		for(int j = 0 ; j < sizeofTago ; j++) {
			int temp = tagoPopulationMatrix[j][0] + tagoPopulationMatrix[j][1] + tagoPopulationMatrix[j][2];
			tagoFitnessValue[j] = temp;
		}
		// �p���g�@��
		for(int j = 0 ; j < sizeofTago ; j++) {
			g1 = tagoPopulationMatrix[j][0];
			g2 = tagoPopulationMatrix[j][1];
			g3 = tagoPopulationMatrix[j][2];
			g4 = tagoPopulationMatrix[j][3];
			g5 = tagoPopulationMatrix[j][4];
			g6 = tagoPopulationMatrix[j][5];
			g7 = tagoPopulationMatrix[j][6];
			g8 = tagoPopulationMatrix[j][7]; 
			judgement = Judge(g1, g2, g3, g4, g5, g6, g7, g8);		
			judgement = judgement * 10000000000.0;
			// fitness value  = �A����ƭ� + �g�@��
			tagoFitnessValue[j] += judgement;			
		}
		// �p��A����ƭȪ���������@ (smaller the better)
		for(int j = 0 ; j < sizeofTago ; j++) {
			// �קK���� (���Ƭ� 0 )
			if(tagoFitnessValue[j] == 0){
				tagoFitnesscPow[j] = 0; 
			}
			else {
				tagoFitnesscPow[j] = (1 / sqrt(tagoFitnessValue[j]) );
			}
		}
		// �p��U���V���� Ef1 �M Ef2
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			for(int k = 0 ; k < sizeofTago ; k++){
				// ef1 �B�z�� tagArr = 0 ��
				if(tagArr[k][j] == 0)
					ef1[j] = ef1[j] + tagoFitnesscPow[k];
				// Ef2 �B�z�� tagArr = 1 ��
				if(tagArr[k][j] == 1)
					ef2[j] = ef2[j] + tagoFitnesscPow[k];
			} 
		}
		// Optimal level
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			if(ef1[j] > ef2[j])
				optimalLevel[j] = 0;
			else
				optimalLevel[j] = 1;
		}
		// �N��]�ȩ�J Optimal Level ����]�ǦC
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			if(optimalLevel[j] == 0)
				tagoFinalGene[j] = tagoGeneArr[0][j];
			else
				tagoFinalGene[j] = tagoGeneArr[1][j];
		} 
		// �N�ҿ� Tago �᪺��]�s�J�ڸs
		gExternalCount++;
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			population[gPopulationCount + gExternalCount - 1][j] = tagoFinalGene[j];
		}
	}
}

void Mutation(double mutationRate, double population[][gNumberOfGenes], int currentIteration, int totalIterations) {
	double randomDouble, lowerLimit, upperLimit;
	int randomInt;
	for(int i = 0 ; i < gPopulationCount + gExternalCount ; i++) {
		for(int j = 0 ; j < gNumberOfGenes ; j++) {
			randomDouble = (double) rand() / (RAND_MAX + 1.0); // ��o 0 ~ 1 �����B�I�ƶü�
			randomInt = rand() % 2; // ��o 0 �� 1 ���H�����
			if(randomDouble < mutationRate) {
				// �i�����
				gExternalCount ++; // ���ܫ���鱵�b�ڸs�᭱(�w�g���ܹL���ڸs�A���A�ѻP����) 
				for(int k = 0 ; k < gNumberOfGenes ; k++){
					population[gPopulationCount + gExternalCount - 1][k] = population[j][k];
				}
				// �V������
				if(randomInt == 0) { 
    				randomDouble = (double) rand() / (RAND_MAX + 1.0);
    				// ����]�d���
					// �]���C�Ӱ�]���d��Ȥ��P�A�G�ݭn�D�Ӱ�]���d���(�W���M�U��)
					lowerLimit = GetLowerLimit(j);
					double temp = population[gPopulationCount + gExternalCount][j];
					// ���ܪ��Z�����H���N���ƼW�[���Y�p
    				double distance = (temp - lowerLimit) * randomDouble * (1 - currentIteration / totalIterations);
    				population[gPopulationCount + gExternalCount][j] -= distance;		
				} 
				else {
					//�V�k����
					randomDouble = (double) rand() / (RAND_MAX + 1.0);
					upperLimit = GetUpperLimit(j);
					double temp = population[gPopulationCount + gExternalCount][j];
    				double distance = (upperLimit - temp) * randomDouble * (1 - currentIteration / totalIterations);
    				population[gPopulationCount + gExternalCount][j] += distance;
    			}
			}
		}
	}
}

void GenerateVariety(double population[][gNumberOfGenes], double fitnessValues[], double fitnessPunishment[]) {
	double g1, g2, g3, g4, g5, g6, g7, g8, judgement;
	int mid = gPopulationCount / 2;
	// �h�˩��H���ͦ�
	// ����: ��e����ӭ��骺�A����ƭȬۦP�ɶi��	
	cout << "�ͦ���]�h�˩�(GenerateVariety)..." << endl;
	if (fitnessValues[mid] == fitnessValues[mid + 1]) {
		for (int i = mid; i < gPopulationCount; i++) {
			// �H���ͦ���]
			g1 = Generate(GetLowerLimit(0), GetUpperLimit(0));
			g2 = Generate(GetLowerLimit(1), GetUpperLimit(1));
			g3 = Generate(GetLowerLimit(2), GetUpperLimit(2));
			g4 = Generate(GetLowerLimit(3), GetUpperLimit(3));
			g5 = Generate(GetLowerLimit(4), GetUpperLimit(4));
			g6 = Generate(GetLowerLimit(5), GetUpperLimit(5));
			g7 = Generate(GetLowerLimit(6), GetUpperLimit(6));
			g8 = Generate(GetLowerLimit(7), GetUpperLimit(7));

			// ��s�ڸs������]�M�A����ƭ�
			population[i][0] = g1;
			population[i][1] = g2;
			population[i][2] = g3;
			population[i][3] = g4;
			population[i][4] = g5;
			population[i][5] = g6;
			population[i][6] = g7;
			population[i][7] = g8;
			judgement = Judge(g1, g2, g3, g4, g5, g6, g7, g8);
			fitnessPunishment[i] = judgement * 10000000000.0;
		}
	} 
}


void VerifyGene(double population[][gNumberOfGenes], double fitnessPunishment[]) {
	double g1, g2, g3, g4, g5, g6, g7, g8, judgement;
	cout << "�ˬd�ڸs��]��(VerifyGene)..." << endl;
	for(int i =  0; i < gPopulationCount + gExternalCount ; i++){
		g1 = population[i][0];
		g2 = population[i][1];
		g3 = population[i][2];
		g4 = population[i][3];
		g5 = population[i][4];
		g6 = population[i][5];
		g7 = population[i][6];
		g8 = population[i][7];
		if(g1 < GetLowerLimit(0) || g1 > GetUpperLimit(0)) {
			g1 = Generate(GetLowerLimit(0), GetUpperLimit(0));
			population[i][0] = g1;
		}
		if(g2 < GetLowerLimit(1) || g2 > GetUpperLimit(1)) {
			g2 = Generate(GetLowerLimit(1), GetUpperLimit(1));
			population[i][1] = g2;
		}
		if(g3 < GetLowerLimit(2) || g3 > GetUpperLimit(2)) {
			g3 = Generate(GetLowerLimit(2), GetUpperLimit(2));
			population[i][2] = g3;
		}
		if(g4 < GetLowerLimit(3) || g4 > GetUpperLimit(3)) {
			g4 = Generate(GetLowerLimit(3), GetUpperLimit(3));
			population[i][3] = g4;
		}
		if(g5 < GetLowerLimit(4) || g5 > GetUpperLimit(4)) {
			g5 = Generate(GetLowerLimit(4), GetUpperLimit(4));
			population[i][4] = g5;
		}
		if(g6 < GetLowerLimit(5) || g6 > GetUpperLimit(5)) {
			g6 = Generate(GetLowerLimit(5), GetUpperLimit(5));
			population[i][5] = g6;
		}
		if(g7 < GetLowerLimit(6) || g7 > GetUpperLimit(6)) {
			g7 = Generate(GetLowerLimit(6), GetUpperLimit(6));
			population[i][6] = g7;
		}
		if(g8 < GetLowerLimit(7) || g8 > GetUpperLimit(7)) {
			g8 = Generate(GetLowerLimit(7), GetUpperLimit(7));
			population[i][7] = g8;
		}
		judgement = Judge(g1, g2, g3, g4, g5, g6, g7, g8);		
		fitnessPunishment[i] = judgement * 1000000000.0;						
	}
}

void init(double fitnessValues[], double& fitnessTotalValues, double fitnessPunishment[], double chosenProbability[], double accumulatedProbabilities[]) {
	memset(fitnessValues, 0, sizeof(double) * gPopulationCount + gExpansionCount);
	fitnessTotalValues = 0;
	memset(fitnessPunishment, 0, sizeof(double) * gPopulationCount + gExpansionCount);
	memset(chosenProbability, 0, sizeof(double) * gPopulationCount);
	memset(accumulatedProbabilities, 0, sizeof(double) * gPopulationCount);
	gExternalCount = 0; // �q�� 2 �����N�}�l�u�O�d gPopulationCount �ӳ̱��񪺵��סC
}

void CompareCurrentFitness(double population[][gNumberOfGenes], double fitnessValues[], int currentIteration, double& targetFitness, double targetGene[], int& targetIteration) {
	if(currentIteration == 0) {
		targetFitness = fitnessValues[0];
		for(int i = 0 ; i < gNumberOfGenes ; i++) {
			targetGene[i] = population[0][i];
		}
		return;
	}
	if(targetFitness > fitnessValues[0]) {
		targetFitness = fitnessValues[0];
		targetIteration = currentIteration;
		for(int i = 0 ; i < gNumberOfGenes ; i++) {
			targetGene[i] = population[0][i];
		}
	}
}

int main() {
	srand(time(NULL));
    double population[gPopulationCount + gExpansionCount][gNumberOfGenes] = {0}; // �ڸs�}�C
	double fitnessValues[gPopulationCount + gExpansionCount], fitnessTotalValues = 0;
	double fitnessPunishment[gPopulationCount + gExpansionCount] = {0}; // �x�s�A���t�ƭȧP�_�᪺��(���T=0�A���~=���j��) 
	double chosenProbability[gPopulationCount] = {0}; // �Q��ܪ����v�C�q�� 2 �����N�}�l�u�O�d gPopulationCount �ӳ̱��񪺵��סC
	double accumulatedProbabilities[gPopulationCount] = {0}; // �ֿn���v
	double mutationRate = 0.1; // ���ܲv 
	int currentIteration = 0; // ��e���N�N��(�ݤj��0)
	const int totalIterations = 30000; // ���N����
	double targetFitness, targetGene[8] = {0};
	int targetIteration = 0;
	
	// �Ĥ@��
    GenerateGene(population, fitnessPunishment); // �ͦ���l��]
	Sortpopulation(population, fitnessValues, fitnessPunishment); // �p���]���A����ƭȡA�ù�ڸs�i��Ƨ�
	CompareCurrentFitness(population, fitnessValues, currentIteration, targetFitness, targetGene, targetIteration); // �����e�A����ƭ� 
	fitnessTotalValues = CalculateFitnessTotal(fitnessValues); // �p���`�A����ƭ�
	CalculateChosenProbability(chosenProbability, fitnessValues, fitnessTotalValues); // �p��Q��ܪ����v
	CalculateAccumulatedProbabilities(accumulatedProbabilities, chosenProbability); // �p��ֿn���v
	Crossover(population, chosenProbability, accumulatedProbabilities); // ��e - ���I���n��e
//	Tago(population, accumulatedProbabilities); // �ϥΥФf���ڸs�i���e
	Mutation(mutationRate, population, currentIteration, totalIterations); // ����
	
    // �ĤG������
	do {
		currentIteration++;
		init(fitnessValues, fitnessTotalValues, fitnessPunishment, chosenProbability, accumulatedProbabilities);
		VerifyGene(population, fitnessPunishment); // �ˬd��]���d��A�íp���g�@��(fitnessPunishment)		
		Sortpopulation(population, fitnessValues, fitnessPunishment);
		CompareCurrentFitness(population, fitnessValues, currentIteration, targetFitness, targetGene, targetIteration);
		cout << "�� " << currentIteration << " �N�̨ξA����ƭ� = " << fitnessValues[0] << endl;
		GenerateVariety(population, fitnessValues, fitnessPunishment); // �ͦ���]�h�˩�
		VerifyGene(population, fitnessPunishment);
		Sortpopulation(population, fitnessValues, fitnessPunishment);
		fitnessTotalValues = CalculateFitnessTotal(fitnessValues); // �p���`�A����ƭ�
		CalculateChosenProbability(chosenProbability, fitnessValues, fitnessTotalValues); // �p��Q��ܪ����v
		CalculateAccumulatedProbabilities(accumulatedProbabilities, chosenProbability); // �p��ֿn���v
		Crossover(population, chosenProbability, accumulatedProbabilities); // ��e - ���I���n��e
		Tago(population, accumulatedProbabilities); // �ϥΥФf���ڸs�i���e
		Mutation(mutationRate, population, currentIteration, totalIterations); // ���� 
	} while(currentIteration < totalIterations);
/*	
	VerifyGene(population, fitnessPunishment);
	Sortpopulation(population, fitnessValues, fitnessPunishment);
	for(int i = 0 ; i < gPopulationCount ; i++) {
		cout << "Gene " << i << " : ";
		for(int j = 0 ; j < gNumberOfGenes ; j++) {
			cout << population[i][j] << "  ";
		}
		cout << endl;
		cout << "�A����� = " << fitnessValues[i] <<endl;
	}
*/	
	cout << "�b�� " << targetIteration << " �N" << endl;
	cout << "��o�̤p���A����� = " << targetFitness << endl;
	cout << "Gene = ";
	for(int i = 0 ; i < gNumberOfGenes ; i++) {
		cout << targetGene[i] << " ";
	}	
	
    return 0;
}

