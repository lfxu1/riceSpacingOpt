/**
 * control morphology of leaf
 * @author weiyang
 *
 */

public class LeafFunction {
	//for leaf, wy131120mark
	public static float parabole(float t)
	{
		float result = ConstPara.a*((float)Math.pow(t, 2)) + ConstPara.b*t;
		return (float)Math.atan(result/t)/ConstPara.rad;
	}

	public static float ellipse(float t)
	{
		float x = ConstPara.xc + ConstPara.d * (float)Math.cos(t*0.2*ConstPara.rad);
		float y = ConstPara.yc + ConstPara.e * (float)Math.sin(t*0.2*ConstPara.rad);
		return (float)Math.atan(y/x)/ConstPara.rad;
	}
	
	//每个叶位的最大叶长
	public static double[] bladeLengthCurve(float maxLength)
	{
		float d = 20;
		int nmax = 13;
		double length[] = new double[16];
		for(int i=0; i<length.length; i++) {
			length[i] = maxLength*Math.pow(Math.E,-(Math.pow((i-nmax),2))/2*d);
		}
		return length;
	}
}
