/**
 * @theory Genetic Algorithm
 * @aim optimize rice structure
 * @author weiyang
 * @date 
 */
import java.io.Serializable;

public class GA implements Serializable {
	
	private Chromosome[] population; //��Ⱥ
	private Chromosome elitist; //��Ӣ����
	private int currentGeneration; //��ǰ�Ż�����

	/**
	 * ���캯��
	 */
	public GA() {
		population = new Chromosome[ConstGaPara.populationSize];
		elitist = new Chromosome();
		currentGeneration = 0;
	}
	
	/**
	 * @Description ��Ⱥ��ʼ�����ڲ���ȡֵ����������������򣬲�����newProMut1��newProMut2
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
	 * @Description ������Ӧ��ֵ����Ⱥ��������
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
	 * @Description ѡ����Ӧ��ֵ���ĸ�����Ϊ��Ӣ����
	 */
	public void elite() {
		try {
			elitist = (Chromosome)population[ConstGaPara.populationSize-1].clone2();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @Description ���̶�ѡ�񷨣���׼����ѡ��
	 * @param population ��Ⱥ
	 * @return ����Ⱥ��ѡ����ĸ����±�
	 */
	public int ratioSelect(Chromosome[] population) {
		int index = 0;
		double sum = 0; //������Ⱥ����Ӧ��
		
		//����������Ⱥ(0~popSize-1)����Ӧ��
		for (int mem = 0; mem < population.length; mem++) {
			sum += population[mem].getFitness();
		}
		//��Ⱥ��ÿ������������Ӧ��ֵ
		for (int mem = 0; mem < population.length; mem++) {
			population[mem].setRfitness(population[mem].getFitness()/sum);
		}
		//��Ⱥ��ÿ��������ۻ���Ӧ��ֵ
		population[0].setCfitness(population[0].getRfitness());
		for (int mem = 1; mem < population.length; mem++) {
			population[mem].setCfitness(population[mem-1].getCfitness() + population[mem].getRfitness());
		}
		population[population.length-1].setCfitness(1);
		//ͨ���ۻ���Ӧ��ֵѡ�����
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
	 * @Description ��������
	 */
	public void evolve() {
		Chromosome[] newpopulation = new Chromosome[ConstGaPara.populationSize]; //�²�������Ⱥ
		Chromosome c1; //�������1
		Chromosome c2; //�������2
		String[] sr = new String[2]; //��Ž��������������¸���Ķ������봮
		
		//�������
		try {
			for(int count=0; count<ConstGaPara.populationSize; count+=2) {
				//ѡ�������������
				c1 = (Chromosome)population[ratioSelect(population)].clone(); //������ѡ�����
				c2 = (Chromosome)population[ratioSelect(population)].clone(); //������ѡ�����
				
				//����λ����
				for(int i=0; i<ConstGaPara.geneNumber; i++) { //һ������һ���������ж��Ƿ���н���
					double x = CodeUtil.random(1000)/1000.0;
					if(x < ConstGaPara.pCro) {
						sr = alienPosCross(c1.getGene()[i], c2.getGene()[i]);//����λ����
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
		
		//�������
		double x;
		for (int i=0; i<ConstGaPara.populationSize; i++) {
			//��һ��
			x = CodeUtil.random(1000)/1000.0;
			if(x<ConstGaPara.newProMut1) {
				for (int j=0; j<ConstGaPara.geneNumber; j++) {
					for(int k=0; k<newpopulation[i].getGene()[j].length(); k++) {
						//�ڶ���
						x = CodeUtil.random(1000)/1000.0;
						if (x < ConstGaPara.newProMut2) {
							char tmp = newpopulation[i].getGene()[j].charAt(k);
							String flag = "1";
							if(tmp == '1') {
								flag = "0";
							}
							//����
							newpopulation[i].setGene(
									newpopulation[i].getGene()[j].substring(0, k) + flag 
									+ newpopulation[i].getGene()[j].substring(k+1), j);
						}
					}
				}
			}
		}
		
		//������Ⱥ
		try {
			for(int i=0; i<ConstGaPara.populationSize; i++){
				population[i] = (Chromosome)newpopulation[i].clone();
			}
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * @Description ��������
	 */
	void compete() {
		//��������
		try {
			if(population[population.length-1].getFitness() < elitist.getFitness()) {
				population[0] = (Chromosome)elitist.clone2();
				sort(); //��Ⱥ����
			} else if(population[population.length-1].getFitness() > elitist.getFitness()) {
				elitist = (Chromosome)population[ConstGaPara.populationSize-1].clone2(); //���¾�Ӣ����
			}
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @Description ����λ����
	 * @param s1 ����1
	 * @param s2 ����2
	 * @return �����õ���������
	 */
	public String[] alienPosCross(String s1, String s2) {
		String[] s = new String[2];
		s[0] = s1;
		s[1] = s2;
		
		int numDif = (int)diversity(s1, s2)*s1.length(); //�����ַ�����ͬλ�ĸ���
		if(numDif <= 0) {
			return s;
		}
		//����ͬλ�ñ��浽croPos[]
		int croPos[][] = new int[numDif][2];
		int index = 0;
		for(int i=0; i<s1.length(); i++) {
			if(s1.charAt(i) != s2.charAt(i)) {
				croPos[index][0] = i; //value
				croPos[index][1] = 0; //flag��0��ʾ����ͬ��1��ʾ��ͬ�����Ѿ����н���
				index++;
			}
		}
		//ѡ�������λ�ý��н���
		if(numDif%2!=0) { //��֤����λ����Ϊż��
			numDif++;
		}
		for(int i=0; i<numDif/2; i++) {
			int point;
			do {
				point = CodeUtil.random(croPos.length); //���ѡ��һ��λ��
			}while(croPos[point][1]== 1);
			int croPoint = croPos[point][0];
			//ֵ����
			String strTmp1 = s[0].substring(0,croPoint) + s[1].charAt(croPoint) + s[0].substring(croPoint+1);
			String strTmp2 = s[1].substring(0,croPoint) + s[0].charAt(croPoint) + s[1].substring(croPoint+1);
			
			s[0] = strTmp1;
			s[1] = strTmp2;
			
			croPos[point][1]= 1; //�ѽ��н���
		}
		return s;
	}

	/**
	 * @Description ������������Ĳ����
	 * @param s1 ����1
	 * @param s2 ����2
	 * @return ������Ĳ����
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