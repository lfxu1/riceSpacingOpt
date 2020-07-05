/**
 * @desc Sigmoid Function of Determinate Growth
 * @author weiyang
 */

public class SigmoidFunction {
	/* wy131125
	 * calculate potential growth rate(pgr) at time t
	 * @para Wmax: the maximal value for weight(w), which is achieved at te
	 * @para te: the time at the end of growth
	 * @para tm: the time at which the maximum growth rate is obtained
	 * @para t: time
	 * @return pgr
	 */
	public static float YinDerivative(float Wmax, float te, float tm, int t)
	{
		//max. growth rate in the linear phase: cm
		float cm = Wmax*((2*te-tm)/(te*(te-tm)))* (float)Math.pow((tm/te), (tm/(te-tm)));
		// the growth rate at time t:
		float pgr = cm*((te-t)/(te-tm))* (float)Math.pow((t/tm), (tm/(te-tm)));
		if (t > te) 
			pgr = 0;
		return pgr;
	}

	/* wy131125
	 * calculate actual growth rate(agr)
	 * @para pgr: caculate from YinDerivative(Wmax, te, tm, t)
	 * @para cpool: cpool of a plant
	 */
	 public static double Growth(float pgr, float sinkDemand, float cpool) {
		float ss_rel = pgr / sinkDemand; 
		//float agr = Math.min(ss_rel * Parameter.INSTANCE.Cpool[id-1], pgr); //wy131121change
		float agr = ss_rel * cpool;//Parameter.INSTANCE.Cpool[id-1]; //wy131124
		agr = Math.min(agr, pgr); //wy131121 
		return agr;//
	 }
	 
	 //caculate agr2
	 public static double Growth2(float pgr, float sinkDemand, float cpool) {
		float ss_rel = pgr / sinkDemand; 
		//float agr = Math.min(ss_rel * Parameter.INSTANCE.Cpool[id-1], pgr); //wy131121change
		float agr2 = ss_rel * cpool; //Parameter.INSTANCE.Cpool[id-1]; //wy131124
		return agr2;
	 }
}
