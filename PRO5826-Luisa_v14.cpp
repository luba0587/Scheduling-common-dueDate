/*******************************************
 * PRO5826 - Trabalho Final - Luísa Cavalcanti
 *
 * Problema de sequenciamento de máquina única:  1 / pj / sum_in_j(WTj*Tj+WEj*Ej) 
 * 	280 Instâncias disponíveis no OR-Library
 * 	data de entrega comum a todas às tarefas, problema restritivo
 *  heurística construtiva
 * 
 *******************************************/

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "randomc.h"
#include "mother.cpp"
#include <chrono>
#include "windows.h"


 
using namespace std;


//classe de medição de tempo
class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const { 
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

//estrutura que guarda o tempo computacional despedido em cada instância
vector < vector< vector<double> > > computationalTimes;

//struct tarefa que guarda o índice do job, seu tempo de processamento e penalidades unitárias de atraso/adiantamento
struct job{
	
	//índice do job
	int index;
	
	//tempo de processamento da tarefa (job)
	double processingTime;
	
	//peso de adiantamento, por unidade de tempo atrasada)
	double earliWeight;
	
	//peso de atraso, por unidade de tempo atrasada)
	double tardiWeight;
	
};

//classe de instâncias do problema
class Instances{
	

	private: 
	
	string inputFileName;	//nome do arquivo texto que contém os dados de entrada
	
	int numProblems;		//número de problemas da instância
	int numJobs;      		//número de jobs de cada problema
	

	public:

	//construtor que cria instância usando o nome do arquivo de entrada
	Instances (string str)
	: inputFileName(str)
	{}
	
	//método que carrega os dados no vetor de jobs de cada problema desta instância (primeiro acessa problema, depois tarefa)
	vector< vector<job> > load();
	
};

//função que inicializa as instâncias com os nomes de arquivo txt
vector<Instances> initializeInstances(){
	
	vector<Instances> all;
	
	Instances tenJobs("sch10.txt");
	all.push_back(tenJobs);
	
	Instances twentyJobs("sch20.txt");
	all.push_back(twentyJobs);
	
	Instances fiftyJobs("sch50.txt");
	all.push_back(fiftyJobs);
	
	Instances hundredJobs("sch100.txt");
	all.push_back(hundredJobs);
	
	Instances twohJobs("sch200.txt");
	all.push_back(twohJobs);
	
	Instances fivehJobs("sch500.txt");
	all.push_back(fivehJobs);
	
	Instances thousJobs("sch1000.txt");
	all.push_back(thousJobs);
	
	return all;
	
}

//função que lê os dados do problema e armazena num vetor de tarefas
vector< vector<job> > Instances::load(){

	//cria vetor de jobs, ainda sem dimensão, que será retornado pela função
	vector< vector<job> > instanceData;
	
	//abre o arquivo de entrada da instância
	ifstream dataFile;
	dataFile.open(inputFileName);
	
	//se der erro na abertura do arquivo, imprime aviso
	if (dataFile.bad()){
		cerr << "Nao foi possivel abrir arquivo de entrada!\n";
	}
	
	//se abrir com sucesso, lê e armazena os dados do arquivo
	else{
		
		//lê o número de problemas da instância
		dataFile >> numProblems;

		//para todo problema...
		for(int p=1; p<= numProblems; p++){
			
			//...cria vetor de job para o problema
			vector<job> problemData;
			
			//...cria tarefa que será armazenada no vetor de jobs do problema
			job read;
			
			//...lê o número de jobs do problema a armazenar
			dataFile >> numJobs;
			
			//...loop que lê os dados conforme numero de tarefas especificada no próprio arquivo (numJobs), 
			//salvando do vetor de job criado no início desse loop (problemData)
			for (int j=1; j<= numJobs; j++){
				
				//indez do job é igual ao contador j do loop
				read.index = j;
				
				//guarda valores de pj, EWj e TWj
				dataFile >> read.processingTime >> read.earliWeight >> read.tardiWeight;
				
				problemData.push_back(read);
				
			}//fim do loop que lê os dados de todas as tarefas do problema
			
			//...guarda dados do problema na matriz (vetor de vetores) instanceData
			instanceData.push_back(problemData);
	
		}//fim do loop que lê os dados de todos os problemas da instância
	
    }//fim do procedimento realizado caso o arquivo de entrada tenha sido aberto com sucesso
    
    dataFile.close();
	
	return instanceData;
	
}

//função que lê os benchmarks do problema e armazena numa estrutura - acesso primeiro inst, depois k, depois h
vector< vector< vector<double> > > loadBenchmarks(vector<Instances> &instancias){
	
	//cria estrutura de dados - primeiro acessa instância, depois problema, depois h
	vector< vector< vector<double> > > benchData;

	//abre o arquivo de entrada da instância
	ifstream dataFile;
	dataFile.open("bench.txt");
	
	//se der erro na abertura do arquivo, imprime aviso
	if (dataFile.bad()){
		cerr << "Nao foi possivel abrir arquivo de entrada!\n";
	}
	
	
	//se abrir com sucesso, lê e armazena os dados do arquivo
	else{

		//variável que armazena o número lido de tarefas da instância
		int n;
		
		//lê uma a uma as instâncias do problema
		for(unsigned inst=0; inst<instancias.size(); inst++){

			//lê o número de problemas da instância
			dataFile >> n;

			//cria estrutura para armazena dados da instância (primeiro acessa problema, depois h)
			vector< vector<double> > instValues;
		
			//para todo problema...
			for(int p=1; p<= 10; p++){
			
				//...cria vetor de job para o problema
				vector<double> problemValues(4);
			
				//...lê valores de benchmark para os 4 valores de h, de 0,2 a 0,8
				dataFile >> problemValues[0] >> problemValues[1] >> problemValues[2]>> problemValues[3];

				//...guarda dados do problema na matriz (vetor de vetores) instanceData
				instValues.push_back(problemValues);
				
			}//fim do loop que lê os dados de benchmark de todos os problemas da instância

			benchData.push_back(instValues);
		
		}//fim do procedimento que percorre todas as instâncias do problema
		

    }//fim do procedimento realizado caso o arquivo de entrada tenha sido aberto com sucesso
   
    dataFile.close();

	return benchData;

}

//função que ordena não-crescentemente jobs de acordo com a razão pj /WEj e caso a razão seja a mesma, ordena não-crescentemente por pj
bool decreasingEarliWeightedProcessing(job lhs, job rhs){
	if (lhs.processingTime/lhs.earliWeight == rhs.processingTime/rhs.earliWeight)
		return lhs.processingTime > rhs.processingTime;
	else 
		return lhs.processingTime/lhs.earliWeight > rhs.processingTime/rhs.earliWeight;
}

//função que ordena não-decrescentemente jobs de acordo com a razão pj /WTj e caso a razão seja a mesma, ordena não-decrescentemente por pj
bool increasingTardiWeightedProcessing(job lhs, job rhs){
	if (lhs.processingTime/lhs.tardiWeight == rhs.processingTime/rhs.tardiWeight)
		return lhs.processingTime < rhs.processingTime;
	else 
		return lhs.processingTime/lhs.tardiWeight < rhs.processingTime/rhs.tardiWeight;
}

//função que ordena não-crescentemente jobs de acordo com WTj e para tarefas com mesmo WTj, ordena por WEj
bool decreasingTardiWeight(job lhs, job rhs){
	if (lhs.tardiWeight == rhs.tardiWeight)
		return lhs.earliWeight < rhs.earliWeight;
	else 
		return lhs.tardiWeight > rhs.tardiWeight;
}

//função que retorna a soma dos tempos de processamento de um vetor de tarefas
double sumProcessTimes(vector<job> &tarefas){
	
	double result=0;
	
	for (unsigned j=0; j<tarefas.size(); j++){
		result += tarefas.at(j).processingTime;	
	}
	
	return result;
	
}

//função que retorna a soma dos pesos de adiantamento de um vetor de tarefas
double sumEarliWeight(vector<job> &tarefas){
	
	double result=0;
	
	for (unsigned j=0; j<tarefas.size(); j++){
		result += tarefas.at(j).earliWeight;
	}
	
	return result;
	
}

//funçaõ que retorna a soma dos pesos de atraso de um vetor de tarefas
double sumTardiWeight(vector<job> &tarefas){
	
	double result=0;
	
	for (unsigned j=0; j<tarefas.size(); j++){
		result += tarefas.at(j).tardiWeight;	
	}
	
	return result;
	
}

//função que roda heurística construtiva
vector<int> runConstructiveHeuristic(double h, vector<job> jobVector, double &t0){
	
	//vetor com índice dos jobs da sequência encontrada
	vector<int> jobOrder;

	//vetor com jobs candidatos a adiantados
	vector<job> earlyCandidates;
	
	//vetor com jobs candidatos a atrasados
	vector<job> lateCandidates;
	
	//calcula e armazena a data de entrega dos jobs do problema e já calcula seu valor
	double dueDate = (int) (h * sumProcessTimes(jobVector) );
	
	//populo os grupos de jobs candidatos a terminar adiantados e atrasados
	for(unsigned j=0; j<jobVector.size(); j++){
		
		if( jobVector.at(j).earliWeight < jobVector.at(j).tardiWeight ) 
			earlyCandidates.push_back( jobVector.at(j) );
		
		else 
			lateCandidates.push_back( jobVector.at(j) );
		
	}
	
	//inicializa variável t
	double t=0;
	
	//verifica o tempo total de processamento dos jobs do grupo E:
	//se todas as tarefas candidatas a terminar adiantadas não couberem antes da data de entrega...
	if (sumProcessTimes(earlyCandidates) > dueDate){
		
		//...ordena o vetor de candidatos decrescentemente pelo peso de atraso
		sort(earlyCandidates.begin(),earlyCandidates.end(),decreasingTardiWeight);
		
		//...tira um a um os elementos com menor peso de atraso até que o tempo de processamento das tarefas remanescentes seja igual ou inferior à due date
		while (sumProcessTimes(earlyCandidates) > dueDate){
			lateCandidates.push_back( earlyCandidates.back() );
			earlyCandidates.pop_back();
		}
			
	}//fim da tratativa para casos em que os elementos de E não cabem antes do prazo
	
	//ordena grupo de candidatos a adiantados não-crescentemente em pj/WEj, ou não-crescentemente em pj quando tal razão for igual
	if (earlyCandidates.size()>0) {
		sort(earlyCandidates.begin(), earlyCandidates.end(), decreasingEarliWeightedProcessing);
	}
	
	
	//ordena grupo de candidatos a atrasados, não-decrescentemente em pj/WTj, (não-decrescentemente em pj como critério de desempate)
	if (lateCandidates.size()>0 ) {
		sort( lateCandidates.begin(), lateCandidates.end(), increasingTardiWeightedProcessing);
	}
	
	//calcula o instante de início de processamento:
	
	//calcula t, que é a soma dos tempos de processamento das tarefas adiantadas
	t = sumProcessTimes(earlyCandidates);
	
	//Se a soma dos pesos de adiantamento for maior que a dos pesos de atraso, não quero adiantar mais as tarefas de E
	if(sumTardiWeight(lateCandidates) < sumEarliWeight(earlyCandidates) )
		t0 = dueDate - t;
	
	//Se fizer sentido adiantar mais as tarefas de E...
	else{
		
		//...se for possível terminar a primeira tarefa de T na data de entrega, o faço
		if ( t + lateCandidates.at(0).processingTime <= dueDate )
			t0 = dueDate - t - lateCandidates.at(0).processingTime;
		
		//...senão, simplesmente começo em t0=0
		else t0=0;
		
	}
	
	//guarda solução encontrada pela heurística construtiva
	for (unsigned j = 0; j < earlyCandidates.size(); j++) {
		jobOrder.push_back( earlyCandidates.at(j).index );
	}
	
	for (unsigned j = 0; j < lateCandidates.size(); j++) {
		jobOrder.push_back( lateCandidates.at(j).index );
	}
	
	return jobOrder;
	
}

//função que calcula o valor de uma solução (sequência e instante de início de processamento) na função objetivo
double calculateFO(vector<int> orderedIndex, vector<job> jobData, double dueDate, double t0){
	
	double resp=0;
	double t=t0;
	
	for(unsigned pos=0; pos<orderedIndex.size(); pos++){
		
		double pj, WEj, WTj;
		
		for(unsigned j=0; j<jobData.size(); j++){
			
			if (jobData.at(j).index == orderedIndex.at(pos)){
				
				pj = jobData.at(j).processingTime;
				
				WEj = jobData.at(j).earliWeight;
				
				WTj = jobData.at(j).tardiWeight;
				
				break;
				
			}
			
		}
		
		t += pj;
		
		if(t <= dueDate) resp += WEj*(dueDate - t);
		
		else resp += WTj*(t-dueDate);
		
	}
	
	return resp;
	
}

//aplica movimento de troca entre uma tarefa adiantada e outra atrasada, dadas as posições e em E e t em T da solução original - retorna 0 se filho é viável e 1 caso contrário
int specificInterchange(double d, vector<job> &sunLate, vector<job> &sunEarly, int e, int t){
	
	//variável de retorno da função - vale 0 se filho é viável e 1 se não foi possível realziar a troca
	int returnValue = 0;

	//calcula o tempo livre antes do prazo tirando a tarefa 'e' de E
	int slackTime = d - sumProcessTimes(sunEarly) + sunEarly.at(e).processingTime;

	//se o tempo disponível antes de d for maior ou igual ao tempo de processamento da tarefa t, faz a troca e registra filho viável
	if (slackTime >= sunLate.at(t).processingTime) {

		//aplica troca entre os dois elementos dessas posições
		job a = sunEarly.at(e);
		sunEarly.at(e) = sunLate.at(t);
		sunLate.at(t) = a;

	}
	

	//caso contrário, executa rotina que atrasa tarefas de E até que a troca seja possível
	else {

		//enquanto a folga disponível antes de d for inferior ao tempo de processamento da tarefa t e houver tarefas em E
		while (slackTime < sunLate.at(t).processingTime && sunEarly.size() > 0 ){

			//verifica se a última tarefa de E é justamente a tarefa e
			if( sunEarly.size()-1 == e){

				//caso positivo, apenas atrasa tarefa e, sem atualização do tempo de folga
				sunLate.push_back(sunEarly.back());
				sunEarly.pop_back();
			
			}

			//caso negativo, atrasa a última tarefa de E e atualiza tempo de folga
			else {

				//insere última tarefa do grupo das adiantadas no grupo das atrasadas e a remove do grupo de adiantadas
				sunLate.push_back(sunEarly.back());
				sunEarly.pop_back();

				//atualiza a folga (desde que a tarefa removida não tenha sido justamente a da posição 'e')
				slackTime += sunLate.back().processingTime;

			}

		}//fim da rotina que remove uma a uma as tarefas de E até que a folga de tempo disponível seja suficiente para encaixar a tarefa t antes da data de entrega

		//se ao fim da rotina não for possível fazer a troca, retorna valor 1
		if (slackTime < sunLate.at(t).processingTime) returnValue = 1;
		
		//caso contrário, faz a troca
		else{
			
			//se a tarefa 'e' ainda estiver no grupo E, faz a troca entre as tarefas 't' e 'e'
			if (e < sunEarly.size()) {
				
				job a = sunEarly.at(e);

				sunEarly.at(e) = sunLate.at(t);
				sunLate.at(t) = a;

			}

			//senão, apenas insere tarefa 't' no grupo de adiantadas
			else {

				sunEarly.push_back(sunLate.at(t));

				sunLate.at(t) = sunLate.back();
				sunLate.pop_back();

			}
		
		}


	}

	return returnValue;
			
}

//aplica movimento que insere uma tarefa especificada por 'posicao' no grupo ao qual ela não pertence originalmente (atrasados ou adiantados) - retorna 0 se o filho é viável, 1 caso contrário
int specificInsertion(double d, vector<job> &sunLate, vector<job> &sunEarly, int posicao){
	
	int returnValue=0;
	
	unsigned n = sunLate.size() + sunEarly.size();
	
	//se a tarefa sorteada é do grupo de adiantadas, a posiciona no grupo das atrasadas
	if ( (unsigned) posicao < sunEarly.size() ){
		
		sunLate.push_back( sunEarly.at(posicao) );
		sunEarly.at(posicao) = sunEarly.back();
		sunEarly.pop_back();
		
	}
	
	//se a tarefa sorteada é do grupo das atrasadas, checa viabilidade do filho e cria nova solução somente se for viável!
	else{
		
		//se a tarefa sorteada não couber no grupo de não-atrasadas, tenta encaixa-la atrasando as tarefas originalmente ñ atrasadas, retorna valor 1 caso ñ seja possível
		if (d - sumProcessTimes(sunEarly) - sunLate.at(n-posicao-1).processingTime < 0){
			
			unsigned maxTrials = sunEarly.size();
			
			//tenta encaixar a tarefa antes do prazo de entrega, atrasando uma a uma as tarefas originalmente no grupo de não-atrasadas com menor peso de atraso
			while(d - sumProcessTimes(sunEarly) - sunLate.at(n-posicao-1).processingTime < 0 && maxTrials>0){
				sunLate.push_back(sunEarly.back());
				sunEarly.pop_back();
				maxTrials -= 1;
			}
			
			//se não conseguir encaixar a tarefa, retorna valor 1
			if (d - sumProcessTimes(sunEarly) - sunLate.at(n-posicao-1).processingTime < 0) 
				returnValue = 1;
			
			//se conseguir, faz a inserção
			else {
				sunEarly.push_back( sunLate.at(n-posicao-1) );
				sunLate.at(n-posicao-1) = sunLate.back();
				sunLate.pop_back();
			}
			
		}
		
		//se a tarefa sorteada couber no grupo de não-atrasadas, faz a inserção
		else{
			sunEarly.push_back( sunLate.at(n-posicao-1) );
			sunLate.at(n-posicao-1) = sunLate.back();
			sunLate.pop_back();
		}
		
		
	}
	
	return returnValue;
	
}


//função que avalia uma solução-filho (candidata a incumbente) dadas as tarefas adiantadas, atrasadas, instante de início de processamento e prazo da instância
double evaluateSun( vector<job> &earlyJobs, vector<job> &lateJobs, double t0, double d){
	
	double result=0;
	double t=t0;
	
	for(unsigned j=0; j< earlyJobs.size(); j++){
		
		t += earlyJobs.at(j).processingTime;
		result += (d-t)*earlyJobs.at(j).earliWeight;
		
	}
	
	for(unsigned j=0; j< lateJobs.size(); j++){
		
		t += lateJobs.at(j).processingTime;
		result += (t-d)*lateJobs.at(j).tardiWeight;
		
	}
	
	return result;
	
}

unsigned contadorInstancias=0;
//função que gera arquivo de saída para acompanhamento da busca local
void generateLocalSearchEvolutionFile ( vector < vector<double> > &dataMatrix){
	
	contadorInstancias += 1;
	
	//nome do arquivo com tabelas de resultados
	ostringstream outputFile;
	outputFile << "localSearch" << contadorInstancias << ".csv";
						
	//cria arquivo de resultados
	ofstream myfile;
	myfile.open (outputFile.str());
		
	//cabeçalho
	myfile << "iteração;valor do filho;melhor valor\n";
	
	//preenche valores das iterações e resultados
	for(unsigned i=0; i<dataMatrix.size(); i++){
		
		for(unsigned j=0; j<dataMatrix.at(i).size(); j++){
		
			myfile << dataMatrix.at(i).at(j) << ";";
				
		}
		
		myfile << "\n";
			
	}
		
	myfile.close();
	
}


//função que realiza a busca local aplicando o movimento de inserção exaustivamente e retorna valor da melhor solução encontrada
double runExaustiveInsertion(vector<int> &sequence, vector<job> jobData, double dueDate, double &t0 ){
	
	//***** parâmetros da busca local
	
	//número máximo de vezes que a busca local é executada, por estagnação ou por qtde de gerações criadas
	unsigned nMax = 200;
	
	//número de filhos de cada geração (tenta trocar todos os jobs de grupo, ou seja lamda=n)
	unsigned lambda = jobData.size(); 

	//***** estrutura de dados usada na busca local
	
	//estrutura que armazena solução pai de uma geração
	vector<job> parentEarlyJobs;
	vector<job> parentLateJobs;
	double parentSolution=0;
	double parent0=t0;
	
	//solução incumbente - guarda o melhor filho já criado entre todas as gerações
	vector<job> incumbentEarlyJobs; 
	vector<job> incumbentLateJobs;
	double incumbentSolution;
	double incumbent0;
	
		
	//estrutura de dados para acompanhar evolução da busca local (só para avaliar evolução da busca local)
	vector < vector<double> >  localSearchEvolution;		//estrutura dados evlução da busca local
	unsigned evaluatedSolutions = 0; 					//contador para o número de soluções avaliadas
	

	//* popula solução pai, já calculando seu valor de FO:
	
	//inicializa t com valor do tempo de processamento do primeiro job da sequência
	unsigned j=0;
	double t= parent0 + jobData.at( sequence.at(j)-1 ).processingTime;
	
	//popula earlyJobs com os jobs adiantados
	while ( t <= dueDate){
		parentEarlyJobs.push_back( jobData.at( sequence.at(j)-1 ) );
		parentSolution += (dueDate-t)*jobData.at( sequence.at(j)-1 ).earliWeight;
		j += 1;
		t += jobData.at( sequence.at(j)-1 ).processingTime;
	}
	
	//popula late jobs com os demais
	while (j < sequence.size() ){
		parentLateJobs.push_back( jobData.at( sequence.at(j)-1 ) );
		parentSolution += (t-dueDate)*jobData.at( sequence.at(j)-1 ).tardiWeight;
		j+=1;
		if (j<sequence.size()) t += jobData.at( sequence.at(j)-1 ).processingTime;
	}//ao fim da rotina, parentSolution possui o valor do pai dessa geração
	
	
	
	
    //***** início do algoritmo da busca local
    
    //variável que guarda o número da geração em que a solução incumbente foi encontrada
	int generationOfBestSolution = 1; 
    
    //solução incumbente recebe solução inicial, que é a primeira solução pai
    incumbentEarlyJobs = parentEarlyJobs;
    incumbentLateJobs = parentLateJobs;
    incumbentSolution = parentSolution;
    incumbent0= parent0;
    
	
    //LOOP QUE CRIA GERAÇÕES DE FILHOS
    for (unsigned n = 1; n<= nMax; n++){
		
		//cria estrutura para armazenar o melhor filho da geração
		vector<job> bestSunEarlyJobs; 
		vector<job> bestSunLateJobs; 
		double bestSunValue=100000000; 
		double bestSunt0=0; 
		
		//LOOP QUE CRIA lambda FILHOS A CADA GERAÇÃO
		for(unsigned sun=1; sun <= lambda; sun++){
			
			//estrutura para acompanhar evolução da busca local
			vector <double> localSearchIterationResult;
			
			//inciliza filho como cópia do pai
			double sunValue=parentSolution;
			double sunt0=0;
			vector<job> sunEarlyJobs = parentEarlyJobs;
			vector<job> sunLateJobs = parentLateJobs;
				
			int viableSun= specificInsertion(dueDate, sunLateJobs, sunEarlyJobs, (sun-1));			
			
			if (viableSun == 0){
					
				//ordena grupo de candidatos a adiantados não-crescentemente em pj/WEj, ou não-crescentemente em pj quando tal razão for igual
				if (sunEarlyJobs.size()>0) {
					sort(sunEarlyJobs.begin(), sunEarlyJobs.end(), decreasingEarliWeightedProcessing);
				}
					
				//ordena grupo de candidatos a atrasados, não-decrescentemente em pj/WTj, (não-decrescentemente em pj como critério de desempate)
				if (sunLateJobs.size()>0 ) {
					sort( sunLateJobs.begin(), sunLateJobs.end(), increasingTardiWeightedProcessing);
				}
					
				//calcula o instante de início de processamento:	
				//calcula t, que é a soma dos tempos de processamento das tarefas adiantadas
				double t = sumProcessTimes(sunEarlyJobs);
					
				//Se a soma dos pesos de adiantamento for maior que a dos pesos de atraso, não quero adiantar mais as tarefas de E
				if(sumTardiWeight(sunLateJobs) < sumEarliWeight(sunEarlyJobs) )
					sunt0 = dueDate - t;
					
				//Se fizer sentido adiantar mais as tarefas de E...
				else{
						
					//...se for possível terminar a primeira tarefa de T na data de entrega, o faço
					if ( t + sunLateJobs.at(0).processingTime <= dueDate )
						sunt0 = dueDate - t - sunLateJobs.at(0).processingTime;
						
					//...senão, simplesmente começo em t0=0					
					else sunt0=0;
						
				}
					
				//***Fim da rotina de ordenação da solução em formato de V e cálculo de t0, igual à construtiva
					
					
				//avalia o valor da função objetivo
				sunValue = evaluateSun( sunEarlyJobs, sunLateJobs, sunt0, dueDate);
					
				//atualiza contador de número de soluções avaliadas
				evaluatedSolutions += 1;
					
				//verifica se o filho é o melhor entre os irmãos já criados
				if (sunValue < bestSunValue){
						
					//update do melhor filho caso positivo
					bestSunValue = sunValue;
						
					bestSunEarlyJobs.clear();
					bestSunLateJobs.clear();
						
					bestSunEarlyJobs = sunEarlyJobs;
					bestSunLateJobs = sunLateJobs;
					bestSunt0 = sunt0;
						
				}
					
				//guarda resultados da iteração (apenas para acompanhamento da evolução da busca local)
				localSearchIterationResult.push_back((double)evaluatedSolutions);
				localSearchIterationResult.push_back(sunValue);
				localSearchIterationResult.push_back(incumbentSolution);
				localSearchEvolution.push_back(localSearchIterationResult);
			
			}//fim da condicional que só leva filhos viáveis em consideração	
			
			
		}//FIM DO LOOP QUE CRIA TODOS OS FILHOS VIÁVEIS A PARTIR DE UM PAI
	
		//verifica se o melhor filho da última geração criada é melhor que a solução incumbente
		if (bestSunValue < incumbentSolution){
				
			//update da melhor solução, caso positivo
			incumbentSolution = bestSunValue;
				
			//atualiza sequência incumbente
			incumbentEarlyJobs.clear();
			incumbentLateJobs.clear();
				
			incumbentEarlyJobs = bestSunEarlyJobs;
			incumbentLateJobs = bestSunLateJobs;
			incumbent0 = bestSunt0;
				
			//grava número de geração da melhor solução
			generationOfBestSolution = n;
				
		}

		//se o melhor filho da última geração não é melhor que o pai dessa geração, finaliza busca local
		else break;
	
		//pai da próxima geração recebe melhor filho dessa geração
		parentEarlyJobs.clear();
		parentLateJobs.clear();
		parentEarlyJobs = bestSunEarlyJobs;
		parentLateJobs = bestSunLateJobs;
		parent0 = bestSunt0;
		
	
	}//FIM DO LOOP QUE CRIA ATÉ nMax GERAÇÕES DE FILHOS
	
	cout << "Melhor geração: " << generationOfBestSolution << " de "<<  nMax << "\n"; 
	
	//gera arquivo com saídas da busca local (apenas para acompanhar evolução)
	generateLocalSearchEvolutionFile(localSearchEvolution);
	
	//guarda solução encontrada pela busca local no vetor sequence
	sequence.clear();
	for (unsigned j = 0; j < incumbentEarlyJobs.size(); j++) {
		sequence.push_back(incumbentEarlyJobs.at(j).index);
	}

	for (unsigned j = 0; j < incumbentLateJobs.size(); j++) {
		sequence.push_back(incumbentLateJobs.at(j).index);
	}
	
	//guarda instante de início de processamento da solução encontrada
	t0 = incumbent0;
	
	//retorna melhor valor de FO encontrado
	return incumbentSolution;
	
}

//função que realiza a busca local aplicando o movimento de troca exaustivamente e retorna valor da melhor solução encontrada
double runExaustiveInterchange(vector<int> &sequence, vector<job> jobData, double dueDate, double &t0 ){
	
	//***** parâmetros da busca local
	
	//número máximo de vezes que a busca local é executada, por estagnação ou por qtde de gerações criadas
	unsigned nMax = 2;

	//***** estrutura de dados usada na busca local
	
	//estrutura que armazena solução pai de uma geração
	vector<job> parentEarlyJobs;
	vector<job> parentLateJobs;
	double parentSolution=0;
	double parent0=t0;
	
	//solução incumbente - guarda o melhor filho já criado entre todas as gerações
	vector<job> incumbentEarlyJobs; 
	vector<job> incumbentLateJobs;
	double incumbentSolution;
	double incumbent0;
	
		
	//estrutura de dados para acompanhar evolução da busca local (só para avaliar evolução da busca local)
	vector < vector<double> >  localSearchEvolution;		//estrutura dados evlução da busca local
	unsigned evaluatedSolutions = 0; 					//contador para o número de soluções avaliadas
	

	//* popula solução pai, já calculando seu valor de FO:
	
	//inicializa t com valor do tempo de processamento do primeiro job da sequência
	unsigned j=0;
	double t= parent0 + jobData.at( sequence.at(j)-1 ).processingTime;
	
	//popula earlyJobs com os jobs adiantados
	while ( t <= dueDate){
		parentEarlyJobs.push_back( jobData.at( sequence.at(j)-1 ) );
		parentSolution += (dueDate-t)*jobData.at( sequence.at(j)-1 ).earliWeight;
		j += 1;
		t += jobData.at( sequence.at(j)-1 ).processingTime;
	}
	
	//popula late jobs com os demais
	while (j < sequence.size() ){
		parentLateJobs.push_back( jobData.at( sequence.at(j)-1 ) );
		parentSolution += (t-dueDate)*jobData.at( sequence.at(j)-1 ).tardiWeight;
		j+=1;
		if (j<sequence.size()) t += jobData.at( sequence.at(j)-1 ).processingTime;
	}//ao fim da rotina, parentSolution possui o valor do pai dessa geração
	
	
    //***** início do algoritmo da busca local
    
    //variável que guarda o número da geração em que a solução incumbente foi encontrada
	int generationOfBestSolution = 1; 
    
    //solução incumbente recebe solução inicial, que é a primeira solução pai
    incumbentEarlyJobs = parentEarlyJobs;
    incumbentLateJobs = parentLateJobs;
    incumbentSolution = parentSolution;
    incumbent0= parent0;
   
	
    //LOOP QUE CRIA GERAÇÕES DE FILHOS
    for (unsigned n = 1; n<= nMax; n++){
		
		//cria estrutura para armazenar o melhor filho da geração
		vector<job> bestSunEarlyJobs; 
		vector<job> bestSunLateJobs; 
		double bestSunValue=100000000; 
		double bestSunt0=0;
		
		//LOOP QUE CRIA E*T FILHOS A CADA GERAÇÃO
		for(unsigned earlyJob=0; earlyJob < parentEarlyJobs.size(); earlyJob++){
			
			//calcula o tempo livre tirando a tarefa 'e' de E
			int slackTime = parent0 - parentEarlyJobs.at(earlyJob).processingTime;
			
			for(unsigned lateJob=0; lateJob < parentLateJobs.size(); lateJob++){
			
				//estrutura para acompanhar evolução da busca local
				vector <double> localSearchIterationResult;
				
				//inciliza filho como cópia do pai
				double sunValue=parentSolution;
				double sunt0=0;
				vector<job> sunEarlyJobs = parentEarlyJobs;
				vector<job> sunLateJobs = parentLateJobs;
				
				//tenta realizar troca entre os elementos earlyJob e lateJob
				int viableSun = specificInterchange(dueDate,sunLateJobs,sunEarlyJobs,earlyJob,lateJob);
				
				
				if (viableSun == 0){
						
					//ordena grupo de candidatos a adiantados não-crescentemente em pj/WEj, ou não-crescentemente em pj quando tal razão for igual
					if (sunEarlyJobs.size()>0) {
						sort(sunEarlyJobs.begin(), sunEarlyJobs.end(), decreasingEarliWeightedProcessing);
					}
						
					//ordena grupo de candidatos a atrasados, não-decrescentemente em pj/WTj, (não-decrescentemente em pj como critério de desempate)
					if (sunLateJobs.size()>0 ) {
						sort( sunLateJobs.begin(), sunLateJobs.end(), increasingTardiWeightedProcessing);
					}
						
					//calcula o instante de início de processamento:	
					//calcula t, que é a soma dos tempos de processamento das tarefas adiantadas
					double t = sumProcessTimes(sunEarlyJobs);
						
					//Se a soma dos pesos de adiantamento for maior que a dos pesos de atraso, não quero adiantar mais as tarefas de E
					if(sumTardiWeight(sunLateJobs) < sumEarliWeight(sunEarlyJobs) )
						sunt0 = dueDate - t;
						
					//Se fizer sentido adiantar mais as tarefas de E...
					else{
							
						//...se for possível terminar a primeira tarefa de T na data de entrega, o faço
						if ( t + sunLateJobs.at(0).processingTime <= dueDate )
							sunt0 = dueDate - t - sunLateJobs.at(0).processingTime;
							
						//...senão, simplesmente começo em t0=0					
						else sunt0=0;
							
					}
						
					//***Fim da rotina de ordenação da solução em formato de V e cálculo de t0, igual à construtiva
						
						
					//avalia o valor da função objetivo
					sunValue = evaluateSun( sunEarlyJobs, sunLateJobs, sunt0, dueDate);
						
					//atualiza contador de número de soluções avaliadas
					evaluatedSolutions += 1;
						
					//verifica se o filho é o melhor entre os irmãos já criados
					if (sunValue < bestSunValue){
							
						//update do melhor filho caso positivo
						bestSunValue = sunValue;
							
						bestSunEarlyJobs.clear();
						bestSunLateJobs.clear();
							
						bestSunEarlyJobs = sunEarlyJobs;
						bestSunLateJobs = sunLateJobs;
						bestSunt0 = sunt0;
							
					}
						
					//guarda resultados da iteração (apenas para acompanhamento da evolução da busca local)
					localSearchIterationResult.push_back((double)evaluatedSolutions);
					localSearchIterationResult.push_back(sunValue);
					localSearchIterationResult.push_back(incumbentSolution);
					localSearchEvolution.push_back(localSearchIterationResult);
				
				}//fim da condicional que só leva filhos viáveis em consideração
				
			}//fim do loop que tenta trocar um job em E por todos os jobs de T
			
		}//fim do loop que tenta trocar todos os jobs de E por todos os de T (GERAÇÃO COMPLETA)
	
	
		//verifica se o melhor filho da última geração criada é melhor que a solução incumbente
		if (bestSunValue < incumbentSolution){
				
			//update da melhor solução, caso positivo
			incumbentSolution = bestSunValue;
				
			//atualiza sequência incumbente
			incumbentEarlyJobs.clear();
			incumbentLateJobs.clear();
				
			incumbentEarlyJobs = bestSunEarlyJobs;
			incumbentLateJobs = bestSunLateJobs;
			incumbent0 = bestSunt0;
				
			//grava número de geração da melhor solução
			generationOfBestSolution = n;
				
		}

		//se o melhor filho da última geração não é melhor que o pai dessa geração, finaliza busca local
		else break;
	
		//pai da próxima geração recebe melhor filho dessa geração
		parentEarlyJobs.clear();
		parentLateJobs.clear();
		parentEarlyJobs = bestSunEarlyJobs;
		parentLateJobs = bestSunLateJobs;
		parent0 = bestSunt0;
		
	
	}//FIM DO LOOP QUE CRIA ATÉ nMax GERAÇÕES DE FILHOS
	
	cout << "Melhor geração: " << generationOfBestSolution << " de "<<  nMax << "\n"; 
	
	//gera arquivo com saídas da busca local (apenas para acompanhar evolução)
	generateLocalSearchEvolutionFile(localSearchEvolution);
	
	//guarda solução encontrada pela busca local no vetor sequence
	sequence.clear();
	for (unsigned j = 0; j < incumbentEarlyJobs.size(); j++) {
		sequence.push_back(incumbentEarlyJobs.at(j).index);
	}

	for (unsigned j = 0; j < incumbentLateJobs.size(); j++) {
		sequence.push_back(incumbentLateJobs.at(j).index);
	}
	
	//guarda instante de início de processamento da solução encontrada
	t0 = incumbent0;
	
	//retorna melhor valor de FO encontrado
	return incumbentSolution;
	
}

//função que cria até nMax gerações compostas de lambda filhos, criados a partir de mutações (inserção e troca) da solução pai. 
//A cada geração, o melhor filho é selecionado como próximo pai, mesmo que haja piora na FO. 
//Se a busca estiver estagnada por stagMax gerações, o pai escolhido é o pior filho da geração.
double runEvolutionaryAlgorithm (vector<int> &sequence, vector<job> jobData,double h, double dueDate, double &t0 ){
	
	//***** parâmetros da busca local
	
	//número máximo de vezes que a busca local é executada, por estagnação ou por qtde de gerações criadas
	unsigned nMax;
	if (jobData.size() <200) nMax = jobData.size();
	else nMax = 100;

	unsigned stagMax=15;
	
	//probabilidade do movimento de inserção
	int insertionProb;
	if ( h < 0.5) insertionProb = 30;
	else insertionProb = 60;
	
	
	//número de filhos de cada geração (em função do tamanho da instância)
	unsigned lambda;
	if (jobData.size() < 100) lambda = jobData.size()*2;
	else {
		if (jobData.size() < 1000) lambda = jobData.size();
		else lambda = jobData.size()/2;
	}


    //semente usada para geração de aleatórios
	//int seed = 5;

	//objetos do gerador de aleatórios
	CRandomMother movementChoice(time(NULL));		//objeto de sorteio usado na escolha dos movimentos
	CRandomMother randPosition(GetTickCount());		//objeto de sorteio usado para selecionar posição aleatória
	
	
	//***** estrutura de dados usada na heurística

	//estrutura que armazena solução pai de uma geração
	vector<job> parentEarlyJobs;
	vector<job> parentLateJobs;
	double parentSolution = 0;
	double parent0 = t0;

	//solução incumbente - guarda o melhor filho já criado entre todas as gerações
	vector<job> incumbentEarlyJobs;
	vector<job> incumbentLateJobs;
	double incumbentSolution;
	double incumbent0;

	//estrutura de dados para acompanhar evolução da busca local (só para avaliar evolução da heurística)
	//vector < vector<double> >  localSearchEvolution;		//estrutura dados evlução auxiliar (só para avaliar evolução)
	
	//contador para o número de soluções avaliadas
	unsigned evaluatedSolutions = 0; 						


	//* popula solução pai, já calculando seu valor de FO:

	//inicializa t com valor do tempo de processamento do primeiro job da sequência
	unsigned j = 0;
	double t = parent0 + jobData.at(sequence.at(j) - 1).processingTime;

	//popula earlyJobs com os jobs não atrasados, na ordem em que são processados
	while (t <= dueDate) {
		parentEarlyJobs.push_back(jobData.at(sequence.at(j) - 1));
		parentSolution += (dueDate - t)*jobData.at(sequence.at(j) - 1).earliWeight;
		j += 1;
		t += jobData.at(sequence.at(j) - 1).processingTime;
	}

	//popula late jobs com os demais, na ordem em que são processados
	while (j < sequence.size()) {
		parentLateJobs.push_back(jobData.at(sequence.at(j) - 1));
		parentSolution += (t - dueDate)*jobData.at(sequence.at(j) - 1).tardiWeight;
		j += 1;
		if (j<sequence.size()) t += jobData.at(sequence.at(j) - 1).processingTime;
	}//ao fim da rotina, parentSolution possui o valor do pai dessa geração




	 //***** início do algoritmo da busca local

	 //variável que guarda o número da geração em que a solução incumbente foi encontrada
	int generationOfBestSolution = 1;

	//variável que conta o número de gerações sem melhoria na solução incumbente
	int stagCounter = 0;


	//solução incumbente recebe solução inicial, que é a primeira solução pai
	incumbentEarlyJobs = parentEarlyJobs;
	incumbentLateJobs = parentLateJobs;
	incumbentSolution = parentSolution;
	incumbent0 = parent0;


	//***** arquivo de saída
	//nome do arquivo com tabelas de resultados
	contadorInstancias++;
	ostringstream outputFile;
	outputFile << "EA-" << contadorInstancias << ".csv";

	//cria arquivo de resultados
	ofstream myfile;
	myfile.open(outputFile.str());

	//cabeçalho
	myfile << "iteracao;geracao;pai;movimento;filho;s*\n";


	//LOOP QUE CRIA nMax GERAÇÕES DE FILHOS
	for (unsigned n = 1; n <= nMax; n++) {

		//cria estrutura para armazenar o melhor filho da geração
		vector<job> bestSunEarlyJobs;
		vector<job> bestSunLateJobs;
		double bestSunValue = 100000000;
		double bestSunt0 = 0;

		//cria estrutura para armazenar o pior filho da geração
		vector<job> worstSunEarlyJobs;
		vector<job> worstSunLateJobs;
		double worstSunValue = 0;
		double worstSunt0 = 0;

		//variáveis que contam o número de vezes que cada movimento foi usado com sucesso
		//int insrtCounter = 0;
		//int chngCounter = 0;

		//LOOP QUE CRIA lambda FILHOS A CADA GERAÇÃO
		for (unsigned sun = 1; sun <= lambda; sun++) {

			//estrutura para acompanhar evolução da busca local
			vector <double> localSearchIterationResult;

			//nome do movimento aplicado
			string movement;

			//inciliza filho como cópia do pai
			double sunValue = parentSolution;
			double sunt0 = 0;
			vector<job> sunEarlyJobs = parentEarlyJobs;
			vector<job> sunLateJobs = parentLateJobs;

			//sorteia o movimento a ser executado na criação do filho
			int move = movementChoice.IRandom(1, 100);
			int viableSun=1;

			//aplica movimento de inserção, se sorteado ou se a solução pai não possui elmentos em E ou em T
			if (move <= insertionProb || sunEarlyJobs.size()<1 || sunLateJobs.size()<1) {
				
				//guarda o nome do movimento aplicado
				movement = "Insercao";

				//escolhe posição do job que terá status alterado
				int p = randPosition.IRandom(0, jobData.size()-1);

				//aplica inserção do job na posição p
				viableSun = specificInsertion(dueDate, sunLateJobs, sunEarlyJobs, p);

			}

			//aplica movimento de troca, quando sorteado e somente se for possível
			else {

				//guarda o nome do movimento aplicado na construção do filho
				movement = "Troca";

				//escolhe duas posições aleatórias em E e em T da solução filho
				int earlyJobPosition = randPosition.IRandom(0, sunEarlyJobs.size() - 1);
				int lateJobPosition = randPosition.IRandom(0, sunLateJobs.size() - 1);

				//tenta realizar troca entre as tarefas sorteadas
				viableSun = specificInterchange(dueDate, sunLateJobs, sunEarlyJobs, earlyJobPosition, lateJobPosition);

			}

			if (viableSun == 0) {

				//ordena grupo de candidatos a adiantados não-crescentemente em pj/WEj, ou não-crescentemente em pj quando tal razão for igual
				if (sunEarlyJobs.size()>0) {
					sort(sunEarlyJobs.begin(), sunEarlyJobs.end(), decreasingEarliWeightedProcessing);
				}

				//ordena grupo de candidatos a atrasados, não-decrescentemente em pj/WTj, (não-decrescentemente em pj como critério de desempate)
				if (sunLateJobs.size()>0) {
					sort(sunLateJobs.begin(), sunLateJobs.end(), increasingTardiWeightedProcessing);
				}

				//calcula o instante de início de processamento:	
				//calcula t, que é a soma dos tempos de processamento das tarefas adiantadas
				double t = sumProcessTimes(sunEarlyJobs);

				//Se a soma dos pesos de adiantamento for maior que a dos pesos de atraso, não quero adiantar mais as tarefas de E
				if (sumTardiWeight(sunLateJobs) < sumEarliWeight(sunEarlyJobs))
					sunt0 = dueDate - t;

				//Se fizer sentido adiantar mais as tarefas de E...
				else {

					//...se for possível terminar a primeira tarefa de T na data de entrega, o faço
					if (t + sunLateJobs.at(0).processingTime <= dueDate)
						sunt0 = dueDate - t - sunLateJobs.at(0).processingTime;

					//...senão, simplesmente começo em t0=0					
					else sunt0 = 0;

				}

				//***Fim da rotina de ordenação da solução em formato de V e cálculo de t0, igual à construtiva


				//avalia o valor da função objetivo
				sunValue = evaluateSun(sunEarlyJobs, sunLateJobs, sunt0, dueDate);

				//atualiza contador de número de soluções avaliadas
				evaluatedSolutions++;

				//verifica se o filho é o melhor entre os irmãos já criados
				if (sunValue < bestSunValue) {

					//update do melhor filho caso positivo
					bestSunValue = sunValue;

					bestSunEarlyJobs.clear();
					bestSunLateJobs.clear();

					bestSunEarlyJobs = sunEarlyJobs;
					bestSunLateJobs = sunLateJobs;
					bestSunt0 = sunt0;

				}

				//verifica se o filho é o pior entre os irmão já avaliados
				if (sunValue > worstSunValue) {

					//update do pior filho caso positivo
					worstSunValue = sunValue;

					worstSunEarlyJobs.clear();
					worstSunLateJobs.clear();

					worstSunEarlyJobs = sunEarlyJobs;
					worstSunLateJobs = sunLateJobs;
					worstSunt0 = sunt0;

				}

				//escreve resultados no arquivo de saída
				myfile << evaluatedSolutions << ";" << n << ";" << parentSolution << ";"
					<< movement << ";" << sunValue << ";" << incumbentSolution << "\n";

				//atualiza contador dos movimentos caso o filho seja melhor que o pai
				//if (sunValue < parentSolution) {

					//if (movement == "Troca") chngCounter++;
					//else insrtCounter++;

				//}


			}//fim da condicional que só leva filhos viáveis em consideração

			//se não tiver sido possível criar filho viável, desconsidera esta iteração
			else sun--;

		}//FIM DO LOOP QUE CRIA lambda FILHOS VIÁVEIS A PARTIR DE UM ÚNICO PAI

		//pai da próxima geração recebe melhor filho dessa geração
		parentSolution = bestSunValue; 
		parentEarlyJobs.clear();
		parentLateJobs.clear();
		parentEarlyJobs = bestSunEarlyJobs;
		parentLateJobs = bestSunLateJobs;
		parent0 = bestSunt0;
		

		//verifica se o melhor filho da última geração criada é melhor que a solução incumbente
		if (bestSunValue < incumbentSolution) {

			//update da melhor solução, caso positivo
			incumbentSolution = bestSunValue;

			//atualiza sequência incumbente
			incumbentEarlyJobs.clear();
			incumbentLateJobs.clear();
			incumbentEarlyJobs = bestSunEarlyJobs;
			incumbentLateJobs = bestSunLateJobs;
			incumbent0 = bestSunt0;

			//grava número de geração da melhor solução
			generationOfBestSolution = n;

			//zera contador de estagnação
			stagCounter = 0;

		}

		//caso contrário, é preciso aumentar o contador de estagnação e verificar se critério de reinício foi atingido
		else {
			
			//incrementa contador de estagnação
			stagCounter++;

			//avalia critério de reinício, caso seja atingido substitui o pai da próx geração pelo pior filho desta última geração
			if (stagCounter >= stagMax) {

				//update do valor da solução pai como pior filho da última geração
				parentSolution = worstSunValue;
				
				//vetores do pai da próxima geração recebe os do pior filho dessa geração
				parentEarlyJobs.clear();
				parentLateJobs.clear();
				parentEarlyJobs = worstSunEarlyJobs;
				parentLateJobs = worstSunLateJobs;
				parent0 = worstSunt0;
				
				//reinicia o contador de estagnação
				stagCounter=0;

			}

		}

		//atualiza a probabilidade de cada movimento
		//insertionProb = insrtCounter / lambda;


	}//FIM DO LOOP QUE CRIA ATÉ nMax GERAÇÕES DE FILHOS

	//fecha arquivo de saída
	myfile.close();

	cout << "Melhor geração: " << generationOfBestSolution << " de " << nMax << "\n";

	//guarda solução encontrada pela busca local no vetor sequence
	sequence.clear();
	for (unsigned j = 0; j < incumbentEarlyJobs.size(); j++) {
		sequence.push_back(incumbentEarlyJobs.at(j).index);
	}

	for (unsigned j = 0; j < incumbentLateJobs.size(); j++) {
		sequence.push_back(incumbentLateJobs.at(j).index);
	}

	//guarda instante de início de processamento da solução encontrada
	t0 = incumbent0;

	//retorna melhor valor de FO encontrado
	return incumbentSolution;
	
}


vector< vector< vector<double> > > calculateAllInstances(vector<Instances> allInstances, vector<double> hFactors){
	
	//estrutura que guarda o valor da FO da melhor solução encontrada pela heurística
	vector< vector< vector<double> > > solutionValues;
	

	//rotina que carrega os dados de cada instância do vetor de instâncias e roda heurística construtiva
	for (unsigned inst=0; inst< allInstances.size(); inst++){
		
		//vetor de vetor de jobs - acesso primeiro o problema, depois a tarefa
		vector< vector<job> > inputData;
		inputData = allInstances.at(inst).load();

		//vetor de vetor de inteiros - índice do job para guardar a sequência obtida pela heurística
		vector< vector< vector<int> > > sequenceOfProblemAndH;
		
		//vetor com tempos de processamento de cada problema para cada valor de h
		vector< vector< double> > compTimeOfProblemAndH;
		
		//vetor de com início de processamento de cada problema para cada h
		vector< vector <double> > startTimeOfProblemAndH;
		
		//vetor de vetor com valores de função objetivo das soluções encontradas(acesso primeiro k, depois h)
		vector< vector<double> > solutionValuesOfInstance;
		
		//para cada problema da instância...
		for (unsigned p=0; p < inputData.size(); p++){
			
			//cria vetor para armazenar o valor de FO do problema, por h adotado
			vector<double> solutionValuesOfProblem;
			
			//vetor com tempo de processamento
			vector <double> compTimeOfH;

			//vetor de vetor de inteiros - índice do job para guardar a sequência obtida pela heurística
			vector< vector<int> > sequenceOfH;
		
			//vetor de com início de processamento de cada problema, por h adotado
			vector <double> startTimeOfH;


			//...e para cada valor de h
			for(unsigned h=0; h< hFactors.size(); h++){
				
				//variável que armazena o tempo de processamento da heurística na instância em questão
				Timer tmr;
			
				//cria vetor para armazenar sequência de jobs encontrada pela heurística 
				vector <int> sequence;
			
				//cria variável para guardar início de processamento da sequência
				double t_start=0;
			
				//chama função da heurística e armazena resultado no vetor sequence e na variável t_start
				sequence = runConstructiveHeuristic(hFactors.at(h), inputData.at(p), t_start);
				
				//computa o resultado na FO da sequência encontrada pela heurística construtiva
				//double temp = calculateFO(sequence, inputData.at(p), (int) (hFactors.at(h)*sumProcessTimes(inputData.at(p))), t_start );
				
				//roda busca local tomando como solução inicial o resultado da construtiva
				double temp = runExaustiveInsertion(sequence, inputData.at(p), (int) (hFactors.at(h)*sumProcessTimes(inputData.at(p))), t_start );
				//double temp = runExaustiveInterchange(sequence,inputData.at(p), (int) (hFactors.at(h)*sumProcessTimes(inputData.at(p))), t_start );
				
				//roda o algoritmo evolucionário usando a saída da construtiva como solução inicial
				//double temp = runEvolutionaryAlgorithm(sequence, inputData.at(p), hFactors.at(h), (int)(hFactors.at(h)*sumProcessTimes(inputData.at(p))), t_start );
				
				//armazena resultados 
				compTimeOfH.push_back(tmr.elapsed());
				sequenceOfH.push_back(sequence);
				startTimeOfH.push_back(t_start);
				solutionValuesOfProblem.push_back(temp);
				
			
			}//fim do loop que percorre todos os valores de h
			
			
			//guardo valores do problema para todos os fatores h nas estruturas apropriadas
			sequenceOfProblemAndH.push_back(sequenceOfH);
			compTimeOfProblemAndH.push_back(compTimeOfH);
			startTimeOfProblemAndH.push_back(startTimeOfH);
			solutionValuesOfInstance.push_back(solutionValuesOfProblem);
	
	
		}//fim do loop que percorre todos os problemas com dado número de jobs

		//guardo os valores da instância
		solutionValues.push_back(solutionValuesOfInstance);
		computationalTimes.push_back(compTimeOfProblemAndH);
		
		
		
		
		//************ROTINA Tabela de Saída***********************************
		
		//nome do arquivo com tabelas de resultados
		ostringstream outputFile;
		outputFile << "sch" << inputData.at(0).size() << ".csv";
						
		//cria arquivo de resultados
		ofstream myfile;
		myfile.open (outputFile.str());
		
		//tabela 2) n, k, h, posição sequência, índice do job
		myfile << "\nInstancia;Problema;Fator h;Posicao;Indice do Job;Start Time\n";
		
		for(unsigned p=0; p < inputData.size(); p++){
			
			for (unsigned h=0; h< hFactors.size(); h++){
				
				for (unsigned j = 0; j < inputData.at(p).size(); j++){
				
					myfile << inst+1 << ";"
						<< p+1 << ";"
						<< hFactors.at(h) << ";"
						<< j+1 << ";"
						<< sequenceOfProblemAndH.at(p).at(h).at(j)  << ";"
						<< startTimeOfProblemAndH.at(p).at(h) << "\n";

				}
				
			}
			
		}
		//************Fim da Tabela de Saída****************************
		
		
		
	
	}//fim do loop que percorre todas as instâncias (todos os arquivos de entrada)

	
	return solutionValues;
	

}



int main(){
	
	//cria vetor com todas as intâncias do problema
	vector< Instances > allInstances = initializeInstances();
	
	//vetor com todos os valores de h utilizados nas instâncias
	vector<double> hFactors;

	//cria estrutura com benchmarks - acesso instância, depois problema, depois h
	vector< vector< vector<double> > > benchmarks;

	//popula vetor hFactors com valores de h
	hFactors.push_back(0.2);
	hFactors.push_back(0.4);
	hFactors.push_back(0.6);
	hFactors.push_back(0.8);
	
	//cria estrutura para guardar os valores obtidos na F.O. - acesso instância, depois problema, depois h
	vector< vector< vector<double> > > solutionValues = calculateAllInstances(allInstances, hFactors);

	//carrega os dados de benchmark
	benchmarks = loadBenchmarks(allInstances);
	
	
	
	//*********************** GERA ARQUIVO DE SAÍDA**********************************
		
	//nome do arquivo com tabelas de resultados
	ostringstream outputFile;
	outputFile << "output-buscaLocal.csv";
						
	//cria arquivo de resultados
	ofstream myfile;
	myfile.open (outputFile.str());
		
	//tabela 1) n,k,h, idle time calculado, valor FO, benchmark, tempo computacional
	myfile << "\nInstancia;Problema;Fator h;ResultadoFO;Benchmark;Tempo Computacional\n";
	
	//preenche valores da tabela 1
	for(unsigned n=0; n < benchmarks.size(); n++){
		
		for(unsigned k=0; k < benchmarks.at(n).size(); k++){
			
			for(unsigned h=0; h < benchmarks.at(n).at(k).size(); h++){
				
				myfile << n+1 << ";"
					<< k+1 << ";"
					<< hFactors.at(h) << ";"
					<< solutionValues.at(n).at(k).at(h) << ";"
					<< benchmarks.at(n).at(k).at(h) << ";"
					<< computationalTimes.at(n).at(k).at(h) << "\n";
				
			}
			
		}
		
	}
		
	myfile.close();
		
	//*********************** fim rotina saída ***********************************/

	//system("pause");

	return 0;
	
}
