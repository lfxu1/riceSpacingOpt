import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * @Description 记录优化过程中的数据处理
 * 
 * methods:
 * clear() 初始化记录文件,清空
 * log(ga) 记录当前遗传算法类ga的染色体数据
 * setTimeLog(fileName, s) 向fileName文件中记录标记为s的当前时间
 * append(fileName, content) 附加信息conten到文件fileName尾
 * overide(fileName, content) 用信息content覆盖文件fileName中的内容
 * div(ga) 计算遗传算法ga当前种群的多样性
 * format(d) 格式化double数d
 * df date数据格式化的格式
 * 
 * @author weiyang
 */
public class ReportUtil {
	
	/**
	 * @Description 优化时，初始化记录文件
	 */
	public static void clear() {
		overide(ConstGaPara.gaLogFile, ""); 
		overide(ConstGaPara.bestPlantFile, ""); 
		overide(ConstGaPara.divFile, "");
		overide(ConstGaPara.timeLog, ""); 
	}
	
	/**
	 * @Description 测试单株时，初始化记录文件
	 */
	 public static void clear2() {
		overide(ConstPara.growthLog, ""); 
	}
		
	/**
	 * @Description 记录优化的过程，包括当前代数、植株编号、优化参数的取值、适应度值
	 * @param ga 遗传算法类对象
	 */
	public static void log(GA ga) {
		append(ConstGaPara.gaLogFile, "G" + ga.getCurrentGeneration() + "------------------------------");
		String report = null;
		for(int i=0; i<ConstGaPara.populationSize; i++) {
			double tmp[] =  new double[ConstGaPara.geneNumber];
			report = "";
			report += "G" + ga.getCurrentGeneration() + "R" + (i+1) + ": ";
			for(int j=0; j<tmp.length; j++) {
				//对第i个参数进行解码，存放到tmp中
				tmp[j] = CodeUtil.decode(ga.getPopulation()[i].getGene()[j], ConstGaPara.lbound[j], ConstGaPara.ubound[j]);
				report += format(tmp[j]) + " ";
			}
			report += " F: " + format2(ga.getPopulation()[i].getFitness());
			//记录优化中每一个个体的数据
			append(ConstGaPara.gaLogFile, report);
			//记录当前代的最优个体
			if(i==ConstGaPara.populationSize-1) {
				append(ConstGaPara.bestPlantFile, report);
			}
		}
	}
	
	/**
	 * @Description 记录时间到日志中
	 * @param fileName 文件名
	 * @param s 标识：start time / finish time
	 */
	public static void setTimeLog(String fileName, String s) {
		String date = df.format(new Date());
		append(fileName, s + ": " + date);
	}
	
	/**
	 * @Description 附加信息到文件尾，不换行，需要换行时传入"\r\n"
	 * @param fileName 文件名
	 * @param content 附加信息
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
	 * @Description 附加信息到文件尾，换行
	 * @param fileName 文件名
	 * @param content 附加信息
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
	 * @Description 用信息覆盖文件中的内容
	 * @param fileName 文件名
	 * @param content 信息
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
	 * @Description 计算种群多样性测度，并将种群多样性测度值添加到文件中
	 * @param ga 优化算法
	 */
	public static void div(GA ga) {
		double div;
		double tmpDiv = 0; //基因位上数字之和
		
		int tmp = 0;
		for(int i=0; i<ConstGaPara.geneNumber; i++) { //第i个基因
			for(int j=0; j<ga.getPopulation()[0].getGene()[i].length(); j++) { //基因上的第j个基因位
				int num1 = 0;
				int num2 = 0;
				for(int k=0; k<ConstGaPara.populationSize; k++) { //第k个染色体
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
	 * @Description 格式化double数
	 * @param d double类型数
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
	
	//date数据格式化
	public static SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

}