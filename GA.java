/**
 * @theory Genetic Algorithm
 * @aim optimize rice structure
 * @author weiyang
 * @date 
 */
import java.io.Serializable;

public class GA implements Serializable {
	
	private Chromosome[] population; //种群
	private Chromosome elitist; //精英个体
	private int currentGeneration; //当前优化代数

	/**
	 * 构造函数
	 */
	public GA() {
		population = new Chromosome[ConstGaPara.populationSize];
		elitist = new Chromosome();
		currentGeneration = 0;
	}
	
	/**
	 * @Description 种群初始化，在参数取值区间内随机产生基因，并计算newProMut1和newProMut2
	 */
	public void initPopulation() {
		//set population parameters
		for(int i=0; i<ConstGaPara.populationSize; i++) {
			population[i] = new Chromosome(ConstGaPara.lbound, ConstGaPara.ubound, ConstGaPara.precision); 
		}
		//set new proMut
		ConstGaPara.newProMut1 = 1-Math.pow((1-ConstGaPara.pMut),population[0].getTotLen());
		ConstGaPara.newProMut2 = ConstGaPara.pMut/ConstGaPara.newProMut1;
	}
	
	/**
	 * @Description 根据适应度值对种群进行排序
	 */
	public void sort() {
		for(int i=1; i<ConstGaPara.populationSize; i++) {
			for(int j=i-1; j>=0; j--) {
				if(population[j].getFitness() > population[j+1].getFitness()) {
					Chromosome tmp = population[j];
					population[j] = population[j+1];
					population[j+1] = tmp;
				}
			}
		}
	}
	
	/**
	 * @Description 选择适应度值最大的个体作为精英个体
	 */
	public void elite() {
		try {
			elitist = (Chromosome)population[ConstGaPara.populationSize-1].clone2();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @Description 轮盘赌选择法，标准比例选择
	 * @param population 种群
	 * @return 从种群中选择出的个体下标
	 */
	public int ratioSelect(Chromosome[] population) {
		int index = 0;
		double sum = 0; //整个种群的适应度
		
		//计算整个种群(0~popSize-1)的适应度
		for (int mem = 0; mem < population.length; mem++) {
			sum += population[mem].getFitness();
		}
		//种群中每个个体的相对适应度值
		for (int mem = 0; mem < population.length; mem++) {
			population[mem].setRfitness(population[mem].getFitness()/sum);
		}
		//种群中每个个体的累积适应度值
		population[0].setCfitness(population[0].getRfitness());
		for (int mem = 1; mem < population.length; mem++) {
			population[mem].setCfitness(population[mem-1].getCfitness() + population[mem].getRfitness());
		}
		population[population.length-1].setCfitness(1);
		//通过累积适应度值选择个体
		double p = CodeUtil.random(1000)/1000.0; //probability to select survivor
		if (p < population[0].getCfitness()) {
			index = 0;
		}else {
			for (int i = 1; i < population.length; i++) {
				if (p >= population[i-1].getCfitness() && p<population[i].getCfitness()) {
					index = i;
				}
			}
		}
		return index;
	}
	
	/**
	 * @Description 进化操作
	 */
	public void evolve() {
		Chromosome[] newpopulation = new Chromosome[ConstGaPara.populationSize]; //新产生的种群
		Chromosome c1; //交叉个体1
		Chromosome c2; //交叉个体2
		String[] sr = new String[2]; //存放交叉操作后的两个新个体的二进制码串
		
		//交叉操作
		try {
			for(int count=0; count<ConstGaPara.populationSize; count+=2) {
				//选择两个交叉个体
				c1 = (Chromosome)population[ratioSelect(population)].clone(); //按比例选择个体
				c2 = (Chromosome)population[ratioSelect(population)].clone(); //按比例选择个体
				
				//相异位交叉
				for(int i=0; i<ConstGaPara.geneNumber; i++) { //一个参数一个参数的判断是否进行交叉
					double x = CodeUtil.random(1000)/1000.0;
					if(x < ConstGaPara.pCro) {
						sr = alienPosCross(c1.getGene()[i], c2.getGene()[i]);//相异位交叉
						c1.setGene(sr[0], i);
						c2.setGene(sr[1], i);
					}
				}
				newpopulation[count] = (Chromosome)c1.clone();
				newpopulation[count+1] = (Chromosome)c2.clone();
			}
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
		//变异操作
		double x;
		for (int i=0; i<ConstGaPara.populationSize; i++) {
			//第一层
			x = CodeUtil.random(1000)/1000.0;
			if(x<ConstGaPara.newProMut1) {
				for (int j=0; j<ConstGaPara.geneNumber; j++) {
					for(int k=0; k<newpopulation[i].getGene()[j].length(); k++) {
						//第二层
						x = CodeUtil.random(1000)/1000.0;
						if (x < ConstGaPara.newProMut2) {
							char tmp = newpopulation[i].getGene()[j].charAt(k);
							String flag = "1";
							if(tmp == '1') {
								flag = "0";
							}
							//变异
							newpopulation[i].setGene(
									newpopulation[i].getGene()[j].substring(0, k) + flag 
									+ newpopulation[i].getGene()[j].substring(k+1), j);
						}
					}
				}
			}
		}
		
		//更新种群
		try {
			for(int i=0; i<ConstGaPara.populationSize; i++){
				population[i] = (Chromosome)newpopulation[i].clone();
			}
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * @Description 竞争操作
	 */
	void compete() {
		//竞争操作
		try {
			if(population[population.length-1].getFitness() < elitist.getFitness()) {
				population[0] = (Chromosome)elitist.clone2();
				sort(); //种群排序
			} else if(population[population.length-1].getFitness() > elitist.getFitness()) {
				elitist = (Chromosome)population[ConstGaPara.populationSize-1].clone2(); //更新精英个体
			}
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @Description 相异位交叉
	 * @param s1 个体1
	 * @param s2 个体2
	 * @return 交叉后得到两个个体
	 */
	public String[] alienPosCross(String s1, String s2) {
		String[] s = new String[2];
		s[0] = s1;
		s[1] = s2;
		
		int numDif = (int)diversity(s1, s2)*s1.length(); //两个字符串不同位的个数
		if(numDif <= 0) {
			return s;
		}
		//将不同位置保存到croPos[]
		int croPos[][] = new int[numDif][2];
		int index = 0;
		for(int i=0; i<s1.length(); i++) {
			if(s1.charAt(i) != s2.charAt(i)) {
				croPos[index][0] = i; //value
				croPos[index][1] = 0; //flag：0表示不相同，1表示相同或者已经进行交叉
				index++;
			}
		}
		//选择相异的位置进行交叉
		if(numDif%2!=0) { //保证相异位个数为偶数
			numDif++;
		}
		for(int i=0; i<numDif/2; i++) {
			int point;
			do {
				point = CodeUtil.random(croPos.length); //随机选择一个位置
			}while(croPos[point][1]== 1);
			int croPoint = croPos[point][0];
			//值交换
			String strTmp1 = s[0].substring(0,croPoint) + s[1].charAt(croPoint) + s[0].substring(croPoint+1);
			String strTmp2 = s[1].substring(0,croPoint) + s[0].charAt(croPoint) + s[1].substring(croPoint+1);
			
			s[0] = strTmp1;
			s[1] = strTmp2;
			
			croPos[point][1]= 1; //已进行交叉
		}
		return s;
	}

	/**
	 * @Description 计算两个个体的差异度
	 * @param s1 个体1
	 * @param s2 个体2
	 * @return 两个体的差异度
	 */
	public double diversity(String s1, String s2) {
		double div = 0;
		if(s1.length() != s2.length()) {
			System.out.println("Two string don't hava the same length in function diversity()!");
			return -1;
		}
		for(int i=0; i<s1.length(); i++) {
			if(s1.charAt(i)==s2.charAt(i)) {
				div += 0;
			} else {
				div += 1;
			}
		}
		div = div/s1.length();
		return div;
	}

	public int getCurrentGeneration() {
		return currentGeneration;
	}

	public void setCurrentGeneration(int generation) {
		this.currentGeneration = generation;
	}

	public Chromosome[] getPopulation() {
		return population;
	}

	public void setPopulation(Chromosome[] population) {
		this.population = population;
	}

	public Chromosome getElitist() {
		return elitist;
	}

	public void setElitist(Chromosome elitist) {
		this.elitist = elitist;
	}
}