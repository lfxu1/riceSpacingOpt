import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * @Description ��¼�Ż������е����ݴ���
 * 
 * methods:
 * clear() ��ʼ����¼�ļ�,���
 * log(ga) ��¼��ǰ�Ŵ��㷨��ga��Ⱦɫ������
 * setTimeLog(fileName, s) ��fileName�ļ��м�¼���Ϊs�ĵ�ǰʱ��
 * append(fileName, content) ������Ϣconten���ļ�fileNameβ
 * overide(fileName, content) ����Ϣcontent�����ļ�fileName�е�����
 * div(ga) �����Ŵ��㷨ga��ǰ��Ⱥ�Ķ�����
 * format(d) ��ʽ��double��d
 * df date���ݸ�ʽ���ĸ�ʽ
 * 
 * @author weiyang
 */
public class ReportUtil {
	
	/**
	 * @Description �Ż�ʱ����ʼ����¼�ļ�
	 */
	public static void clear() {
		overide(ConstGaPara.gaLogFile, ""); 
		overide(ConstGaPara.bestPlantFile, ""); 
		overide(ConstGaPara.divFile, "");
		overide(ConstGaPara.timeLog, ""); 
	}
	
	/**
	 * @Description ���Ե���ʱ����ʼ����¼�ļ�
	 */
	 public static void clear2() {
		overide(ConstPara.growthLog, ""); 
	}
		
	/**
	 * @Description ��¼�Ż��Ĺ��̣�������ǰ������ֲ���š��Ż�������ȡֵ����Ӧ��ֵ
	 * @param ga �Ŵ��㷨�����
	 */
	public static void log(GA ga) {
		append(ConstGaPara.gaLogFile, "G" + ga.getCurrentGeneration() + "------------------------------");
		String report = null;
		for(int i=0; i<ConstGaPara.populationSize; i++) {
			double tmp[] =  new double[ConstGaPara.geneNumber];
			report = "";
			report += "G" + ga.getCurrentGeneration() + "R" + (i+1) + ": ";
			for(int j=0; j<tmp.length; j++) {
				//�Ե�i���������н��룬��ŵ�tmp��
				tmp[j] = CodeUtil.decode(ga.getPopulation()[i].getGene()[j], ConstGaPara.lbound[j], ConstGaPara.ubound[j]);
				report += format(tmp[j]) + " ";
			}
			report += " F: " + format2(ga.getPopulation()[i].getFitness());
			//��¼�Ż���ÿһ�����������
			append(ConstGaPara.gaLogFile, report);
			//��¼��ǰ�������Ÿ���
			if(i==ConstGaPara.populationSize-1) {
				append(ConstGaPara.bestPlantFile, report);
			}
		}
	}
	
	/**
	 * @Description ��¼ʱ�䵽��־��
	 * @param fileName �ļ���
	 * @param s ��ʶ��start time / finish time
	 */
	public static void setTimeLog(String fileName, String s) {
		String date = df.format(new Date());
		append(fileName, s + ": " + date);
	}
	
	/**
	 * @Description ������Ϣ���ļ�β�������У���Ҫ����ʱ����"\r\n"
	 * @param fileName �ļ���
	 * @param content ������Ϣ
	 */
	public static void append2(String fileName, String content) {
		try {
		    FileWriter writer = new FileWriter(fileName, true);
		    writer.write(content);
		    writer.close();
		} catch (IOException e) {
		    e.printStackTrace();
		}
	}
	
	/**
	 * @Description ������Ϣ���ļ�β������
	 * @param fileName �ļ���
	 * @param content ������Ϣ
	 */
	public static void append(String fileName, String content) {
		try {
		    FileWriter writer = new FileWriter(fileName, true);
		    writer.write(content + "\r\n");
		    writer.close();
		} catch (IOException e) {
		    e.printStackTrace();
		}
	}
    
	/**
	 * @Description ����Ϣ�����ļ��е�����
	 * @param fileName �ļ���
	 * @param content ��Ϣ
	 */
	public static void overide(String fileName, String content) {
		try {
		    FileWriter writer = new FileWriter(fileName);
		    writer.write(content);
		    writer.close();
		} catch (IOException e) {
		    e.printStackTrace();
		}
	}
    
	/** 
	 * @Description ������Ⱥ�����Բ�ȣ�������Ⱥ�����Բ��ֵ��ӵ��ļ���
	 * @param ga �Ż��㷨
	 */
	public static void div(GA ga) {
		double div;
		double tmpDiv = 0; //����λ������֮��
		
		int tmp = 0;
		for(int i=0; i<ConstGaPara.geneNumber; i++) { //��i������
			for(int j=0; j<ga.getPopulation()[0].getGene()[i].length(); j++) { //�����ϵĵ�j������λ
				int num1 = 0;
				int num2 = 0;
				for(int k=0; k<ConstGaPara.populationSize; k++) { //��k��Ⱦɫ��
					tmp = Integer.parseInt(""+ga.getPopulation()[k].getGene()[i].charAt(j));
					num1 += tmp;
					num2 += (1-tmp);
				}
				tmpDiv += (Math.max(num1,num2)-Math.min(num1,num2));
				//append(ga.getDivFile(),"Math.max(num1,num2):" + Math.max(num1,num2) + " Math.min(num1,num2):" + Math.min(num1,num2));
				//append(ga.getDivFile(),"tmpDiv" + tmpDiv);
			}
		}
		//append(ga.getDivFile(),"tmpDiv:" + tmpDiv + " ln:" + ga.getPopulation()[0].getTotLen() * ga.getPopSize());
		div = 1-tmpDiv/(ga.getPopulation()[0].getTotLen() * ConstGaPara.populationSize);
		append(ConstGaPara.divFile, "G" + ga.getCurrentGeneration() + ",div: " + format2(div));
		tmpDiv = 0;
	}
	
	/**
	 * @Description ��ʽ��double��
	 * @param d double������
	 */
	public static String format(double d) {
		DecimalFormat df;
		df= new DecimalFormat("###0");
		return df.format(d);
	}
	
	public static String format2(double d) {
		DecimalFormat df;
		df= new DecimalFormat("###0.00");
		return df.format(d);
	}
	
	//date���ݸ�ʽ��
	public static SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

}