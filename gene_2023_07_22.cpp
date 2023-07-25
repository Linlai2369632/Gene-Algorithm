#include <iostream>
#include <time.h>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

const int gPopulationCount = 80; // 族群數量
const int gNumberOfGenes = 8;
const int gExpansionCount = 10000; // 擴充數量。此變數應只出現在其他變數宣告裡。
int gExternalCount = 0; // 交叉和突變時紀錄額外增加的數量

double Generate(double min, double max) {
    // 生成 min 到 max 之間的隨機數
    double c;
    double Beta;
    Beta = (double) rand() / (RAND_MAX + 1.0);
    c = min + Beta * (max - min);
    c = round(c * 1000000) / 1000000.0; // 只保留前6位，後面捨棄。
    return c;
}

double GetUpperLimit(int gene) {
	// 基因上限
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
	// 基因下限
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
    // 計算判斷式錯誤總距離
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

        // 計算判斷式錯誤總距離 
        judgement = Judge(g1, g2, g3, g4, g5, g6, g7, g8);
        fitnessPunishment[i] = judgement * 10000000000.0; // 若有誤差，將其誤差值放大

        // 存入陣列 
        population[i][0] = g1; 
        population[i][1] = g2; 
        population[i][2] = g3; 
        population[i][3] = g4; 
        population[i][4] = g5; 
        population[i][5] = g6; 
        population[i][6] = g7; 
        population[i][7] = g8;
    }
	// 輸出目前族群 
	cout << "初代基因: \n";
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
	double sortedFitnessValues[gPopulationCount + gExpansionCount];    // 儲存排序後的陣列
	int sortedIndices[gPopulationCount + gExpansionCount] = {0};       // 儲存排序後所對應陣列位置
	
	for(int i = 0; i < gPopulationCount + gExternalCount ; i++){
		fitnessValues[i] = population[i][0] +  population[i][1] + population[i][2] + fitnessPunishment[i];  
	}
	
	memcpy(sortedFitnessValues, fitnessValues, sizeof(double) * totalPopulationSize);
	sort(sortedFitnessValues, sortedFitnessValues + totalPopulationSize);
	// 記錄所對應的族群位置
	for(int i = 0 ; i < totalPopulationSize ; i++){
		for(int j = 0; j < gPopulationCount + gExternalCount ; j++){
			if(sortedFitnessValues[i] == fitnessValues[j]){
				sortedIndices[i] = j;
			}
		}
	} 
	//將族群排序後的位置更新
	double populationTemp[gPopulationCount + gExpansionCount][gNumberOfGenes] = {0}; // 暫存排序後族群陣列
	for(int i = 0 ; i < totalPopulationSize ; i++){
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			int temp = sortedIndices[i];
			populationTemp[i][j] = population[temp][j];
		}
	}
	//將位置更新後的暫存陣列 重新塞回初始族群陣列 
	for(int i = 0 ; i < totalPopulationSize ; i++){
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			population[i][j] = populationTemp[i][j];
		}
	} 
	//將排序後的適應函數值，重新放回初始適應函數值陣列 
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
            selectedForCrossover = gPopulationCount - selectedForCrossover + 1; // 選小塊的
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
    // 單點內積函數
    double ret;
    ret = B * x1 + (1 - B) * x2;
    ret = (int)(ret * 1000000) / 1000000.0;
    return ret;
}

double singlePointInnerProductY(double B, double y1, double y2) {
    // 單點內積函數 (Y值)
    double ret;
    ret = (1 - B) * y1 + B * y2;
    ret = (int)(ret * 1000000) / 1000000.0;
    return ret;
}

void Crossover(double population[][gNumberOfGenes], double chosenProbability[], double accumulatedProbabilities[]) {
	int crossRate = static_cast<int>(gPopulationCount * 0.8);
	int selectedForCrossover1, selectedForCrossover2; // 被選擇做交叉的個體
	int selectedCrossGene; // 被選擇做交叉的基因
	double gene1, gene2; // 交叉的基因
	for(int i = 0 ; i < crossRate ; i++) {
		// 比較累積機率，選出做交叉的個體
		selectedForCrossover1 = SelectIndividualForCrossover(accumulatedProbabilities);
		selectedForCrossover2 = SelectIndividualForCrossover(accumulatedProbabilities);
    	
    	selectedCrossGene = rand() % gNumberOfGenes;
    	gene1 = population[selectedForCrossover1][selectedCrossGene];
    	gene2 = population[selectedForCrossover2][selectedCrossGene];
		
		// 要交叉的個體接在族群尾端
        AddIndividualToPopulation(population, selectedForCrossover1);
        AddIndividualToPopulation(population, selectedForCrossover2);
		
		double Beta = static_cast<double>(rand()) / (RAND_MAX + 1.0); // 取的 0 ~ 1 之間的浮點數
		// 內積(單點內積)  (100+2-2=100)
		population[gNumberOfGenes + gExternalCount - 2][selectedCrossGene] = singlePointInnerProductX(Beta, gene1, gene2);
		// 內積(單點內積)  (100+2-1=101)
    	population[gNumberOfGenes + gExternalCount - 1][selectedCrossGene] = singlePointInnerProductY(Beta, gene1, gene2);
		// 單點後所有元素互換
		for (int j = selectedCrossGene; j < gNumberOfGenes; j++) {
			double temp = population[selectedForCrossover1][j];
			population[selectedForCrossover1][j] = population[selectedForCrossover2][j];
			population[selectedForCrossover2][j] = temp;
		}
	}
}

void Tago(double population[][gNumberOfGenes], double accumulatedProbabilities[]) {
	// 此田口表只適用於基因數目等於 8 時
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
	int selectedForCrossover1, selectedForCrossover2; // 被選擇做交叉的個體
	double g1, g2, g3, g4, g5, g6, g7, g8, judgement;
	double tagoGeneArr[2][gNumberOfGenes]; // 田口輪盤
	double tagoPopulationMatrix[sizeofTago][gNumberOfGenes] = {0}; // 田口陣列 
	double tagoFitnessValue[sizeofTago] = {0}; // 田口適應函數值
	double tagoFitnesscPow[sizeofTago] = {0}; // 田口 FitnessValue 平方分之一 
	double ef1[gNumberOfGenes] = {0};
	double ef2[gNumberOfGenes] = {0};
	int optimalLevel[gNumberOfGenes] = {0};
	double tagoFinalGene[gNumberOfGenes] = {0}; // 田口最後基因值
	for(int i = 0 ; i < crossRate ; i++) {
		// 比較累積機率，選出做交叉的個體
		selectedForCrossover1 = SelectIndividualForCrossover(accumulatedProbabilities);
		selectedForCrossover2 = SelectIndividualForCrossover(accumulatedProbabilities);
		// 所選 2 組存入陣列
		for(int j = 0 ; j < gNumberOfGenes ; j++) {
        	tagoGeneArr[0][j] = population[selectedForCrossover1][j];
			tagoGeneArr[1][j] = population[selectedForCrossover2][j];
        }
		//將基因值帶入田口表
    	for(int j = 0 ; j < sizeofTago ; j++) {
    		for(int k = 0 ; k < gNumberOfGenes ; k++){
    			int temp = tagArr[j][k]; 
    			tagoPopulationMatrix[j][k] = tagoGeneArr[temp][k];
			}
		}
		// 計算 fitness value (適應函數值)
		for(int j = 0 ; j < sizeofTago ; j++) {
			int temp = tagoPopulationMatrix[j][0] + tagoPopulationMatrix[j][1] + tagoPopulationMatrix[j][2];
			tagoFitnessValue[j] = temp;
		}
		// 計算懲罰值
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
			// fitness value  = 適應函數值 + 懲罰值
			tagoFitnessValue[j] += judgement;			
		}
		// 計算適應函數值的平方分之一 (smaller the better)
		for(int j = 0 ; j < sizeofTago ; j++) {
			// 避免溢位 (除數為 0 )
			if(tagoFitnessValue[j] == 0){
				tagoFitnesscPow[j] = 0; 
			}
			else {
				tagoFitnesscPow[j] = (1 / sqrt(tagoFitnessValue[j]) );
			}
		}
		// 計算各條染色體 Ef1 和 Ef2
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			for(int k = 0 ; k < sizeofTago ; k++){
				// ef1 處理當 tagArr = 0 時
				if(tagArr[k][j] == 0)
					ef1[j] = ef1[j] + tagoFitnesscPow[k];
				// Ef2 處理當 tagArr = 1 時
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
		// 將基因值放入 Optimal Level 的基因序列
		for(int j = 0 ; j < gNumberOfGenes ; j++){
			if(optimalLevel[j] == 0)
				tagoFinalGene[j] = tagoGeneArr[0][j];
			else
				tagoFinalGene[j] = tagoGeneArr[1][j];
		} 
		// 將所選 Tago 後的基因存入族群
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
			randomDouble = (double) rand() / (RAND_MAX + 1.0); // 獲得 0 ~ 1 之間浮點數亂數
			randomInt = rand() % 2; // 獲得 0 或 1 的隨機整數
			if(randomDouble < mutationRate) {
				// 進行突變
				gExternalCount ++; // 突變後個體接在族群後面(已經突變過的族群，不再參與突變) 
				for(int k = 0 ; k < gNumberOfGenes ; k++){
					population[gPopulationCount + gExternalCount - 1][k] = population[j][k];
				}
				// 向左突變
				if(randomInt == 0) { 
    				randomDouble = (double) rand() / (RAND_MAX + 1.0);
    				// 取基因範圍值
					// 因為每個基因的範圍值不同，故需要求該基因的範圍值(上限和下限)
					lowerLimit = GetLowerLimit(j);
					double temp = population[gPopulationCount + gExternalCount][j];
					// 突變的距離應隨迭代次數增加而縮小
    				double distance = (temp - lowerLimit) * randomDouble * (1 - currentIteration / totalIterations);
    				population[gPopulationCount + gExternalCount][j] -= distance;		
				} 
				else {
					//向右突變
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
	// 多樣性隨機生成
	// 條件: 當前面兩個個體的適應函數值相同時進行	
	cout << "生成基因多樣性(GenerateVariety)..." << endl;
	if (fitnessValues[mid] == fitnessValues[mid + 1]) {
		for (int i = mid; i < gPopulationCount; i++) {
			// 隨機生成基因
			g1 = Generate(GetLowerLimit(0), GetUpperLimit(0));
			g2 = Generate(GetLowerLimit(1), GetUpperLimit(1));
			g3 = Generate(GetLowerLimit(2), GetUpperLimit(2));
			g4 = Generate(GetLowerLimit(3), GetUpperLimit(3));
			g5 = Generate(GetLowerLimit(4), GetUpperLimit(4));
			g6 = Generate(GetLowerLimit(5), GetUpperLimit(5));
			g7 = Generate(GetLowerLimit(6), GetUpperLimit(6));
			g8 = Generate(GetLowerLimit(7), GetUpperLimit(7));

			// 更新族群中的基因和適應函數值
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
	cout << "檢查族群基因中(VerifyGene)..." << endl;
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
	gExternalCount = 0; // 從第 2 次迭代開始只保留 gPopulationCount 個最接近的答案。
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
    double population[gPopulationCount + gExpansionCount][gNumberOfGenes] = {0}; // 族群陣列
	double fitnessValues[gPopulationCount + gExpansionCount], fitnessTotalValues = 0;
	double fitnessPunishment[gPopulationCount + gExpansionCount] = {0}; // 儲存適應含數值判斷後的值(正確=0，錯誤=極大值) 
	double chosenProbability[gPopulationCount] = {0}; // 被選擇的機率。從第 2 次迭代開始只保留 gPopulationCount 個最接近的答案。
	double accumulatedProbabilities[gPopulationCount] = {0}; // 累積機率
	double mutationRate = 0.1; // 突變率 
	int currentIteration = 0; // 當前迭代代數(需大於0)
	const int totalIterations = 30000; // 迭代次數
	double targetFitness, targetGene[8] = {0};
	int targetIteration = 0;
	
	// 第一輪
    GenerateGene(population, fitnessPunishment); // 生成初始基因
	Sortpopulation(population, fitnessValues, fitnessPunishment); // 計算基因的適應函數值，並對族群進行排序
	CompareCurrentFitness(population, fitnessValues, currentIteration, targetFitness, targetGene, targetIteration); // 比較當前適應函數值 
	fitnessTotalValues = CalculateFitnessTotal(fitnessValues); // 計算總適應函數值
	CalculateChosenProbability(chosenProbability, fitnessValues, fitnessTotalValues); // 計算被選擇的機率
	CalculateAccumulatedProbabilities(accumulatedProbabilities, chosenProbability); // 計算累積機率
	Crossover(population, chosenProbability, accumulatedProbabilities); // 交叉 - 單點內積交叉
//	Tago(population, accumulatedProbabilities); // 使用田口表對族群進行交叉
	Mutation(mutationRate, population, currentIteration, totalIterations); // 突變
	
    // 第二輪之後
	do {
		currentIteration++;
		init(fitnessValues, fitnessTotalValues, fitnessPunishment, chosenProbability, accumulatedProbabilities);
		VerifyGene(population, fitnessPunishment); // 檢查基因的範圍，並計算懲罰值(fitnessPunishment)		
		Sortpopulation(population, fitnessValues, fitnessPunishment);
		CompareCurrentFitness(population, fitnessValues, currentIteration, targetFitness, targetGene, targetIteration);
		cout << "第 " << currentIteration << " 代最佳適應函數值 = " << fitnessValues[0] << endl;
		GenerateVariety(population, fitnessValues, fitnessPunishment); // 生成基因多樣性
		VerifyGene(population, fitnessPunishment);
		Sortpopulation(population, fitnessValues, fitnessPunishment);
		fitnessTotalValues = CalculateFitnessTotal(fitnessValues); // 計算總適應函數值
		CalculateChosenProbability(chosenProbability, fitnessValues, fitnessTotalValues); // 計算被選擇的機率
		CalculateAccumulatedProbabilities(accumulatedProbabilities, chosenProbability); // 計算累積機率
		Crossover(population, chosenProbability, accumulatedProbabilities); // 交叉 - 單點內積交叉
		Tago(population, accumulatedProbabilities); // 使用田口表對族群進行交叉
		Mutation(mutationRate, population, currentIteration, totalIterations); // 突變 
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
		cout << "適應函數 = " << fitnessValues[i] <<endl;
	}
*/	
	cout << "在第 " << targetIteration << " 代" << endl;
	cout << "獲得最小的適應函數 = " << targetFitness << endl;
	cout << "Gene = ";
	for(int i = 0 ; i < gNumberOfGenes ; i++) {
		cout << targetGene[i] << " ";
	}	
	
    return 0;
}

