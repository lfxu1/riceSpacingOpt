/**
 * @descr dynamic sun with changement of position and radiation
 * @author weiyang 20131127
 */
 
public class SunFunction {
	
	/**calculate the position and radiation change of sun
	 * position: An algorithm for the computation of the solar position,Roberto Grena,2007
	 * radiation: 太阳辐射计算讲座大全（王炳忠1-5和解答）,1999
	 * @param tt: time
	 * @param dd: day
	 * @param mm: month
	 * @param yy: year
	 * @param delta_t: difference between universal time and terrestrial time 世界时和地球时差异
	 * @param observerLatitude: 纬度
	 * @param observerLongitude: 经度
	 * @param pressure：气压
	 * @param temperature: 温度 period is used for the index of ConstPara.TA[]
	 * @return position
	 * @return radiation
	 * @auther weiyang
	 */
	/*
	public static float[] calcSun (double tt, int dd, int mm, int yy, 
					 double delta_t, double observerLatitude, double observerLongitude, 
					 double pressure, double temperature)*/
	public static float[] getSun( int yy, int mm, int dd, int hh)
					 //double pressure, double temperature)
	{ 
		double jd;
		double jde;
		double heliocLongitude;
		double geocSolarLongitude;
		double rightAscension;
		double declination;
		
		double hourAngle;
		double topocRightAscension;
		double topocDeclination;
		double topocHourAngle;
		double elevationNoRefrac;
		double refractionCorrection;
		double altitude; //高度角： 入射光与地面夹角
		double zenith; //天顶角：入射光与地面法线夹角
		double azimuth; //方位角：从某点的指北方向线起，依顺时针方向到目标方向线之间的水平夹角。
		
		double ut;
		if(hh>=1||hh<=8) {
			ut = hh - 8 +24;
			dd --;
		} else {
			ut = hh - 8;
		}
		//calculate of JD and JDE
		double dYear, dMonth;
		if(mm<=2) {
			dYear = (double)yy - 1.0;
			dMonth = (double)mm + 12.0;
		}
		else {
			dYear = (double)yy;
			dMonth = (double)mm;
		}
		
		double jd_t = (double)Math.floor(365.25*(dYear-2000))
			+(double)Math.floor(30.6001*(dMonth+1))
			+(double)dd + ut/24.0 - 1158.5;
		
		double t = jd_t + ConstPara.delta_t/86400;
		
		//period is used for the index of ConstPara.TA[]
		int startDayNum = TimeFunction.whichDayOfTheYear(ConstPara.YEAR, ConstPara.MONTH, ConstPara.DAY);
		int currentDayNum = TimeFunction.whichDayOfTheYear(yy, mm, dd);
		int period = currentDayNum - startDayNum + 1;
		
		//standard JD and JDE (useless for the computation, they are computed for completeness)
		jde = t + 2452640;
		jd = jd_t + 2452640;
		 
		//heliocentric longitude
		//linear increase + annual harmonic
		double ang = 1.72019e-2*t - 0.0563;
		heliocLongitude = 1.740940 + 1.7202768683e-2*t + 3.34118e-2*Math.sin(ang) + 3.488e-4*Math.sin(2*ang);
		//moon perturbation
		heliocLongitude += 3.13e-5*Math.sin(2.127730e-1*t-0.585);
		//harmonic correction
		heliocLongitude += 1.26e-5*Math.sin(4.243e-3*t+1.46)
			+ 2.35e-5*Math.sin(1.0727e-2*t+0.72)
			+ 2.76e-5*Math.sin(1.5799e-2*t+2.35)
			+ 2.75e-5*Math.sin(2.1551e-2*t-1.98)
			+ 1.26e-5*Math.sin(3.1490e-2*t-0.80);
		//plynoimal correction
		double t2 = 0.001*t;
		heliocLongitude += ((( -2.30796e-7*t2 + 3.7976e-6)*t2 - 2.0458e-5)*t2 + 3.976e-5)*t2*t2;
		//end heliocentric longitude calculation
		
		//correctin to geocentric longitude due to nutation
		double delta_psi = 8.33e-5*Math.sin(9.252e-4*t - 1.173);
		//earth axis inclination
		double epsilon = -6.21e-9*t + 0.409086 + 4.46e-5*Math.sin(9.252e-4*t+0.397);
		//geocentric global solar coordinates:
		//geocentric solar longitude
		geocSolarLongitude = heliocLongitude + Math.PI + delta_psi - 9.932e-5; 
		//geocentric right ascension
		//double sLambda = Math.sin(geocSolarLongitude);
		rightAscension = Math.atan2(Math.sin(geocSolarLongitude)*Math.cos(epsilon), Math.cos(geocSolarLongitude));
		
		//declination
		declination = Math.asin(Math.sin(epsilon)*Math.sin(geocSolarLongitude));
		//local hour angle of the sun
		hourAngle = 6.30038809903*jd_t + 4.8824623 + delta_psi*0.9174 + 
								ConstPara.observerLongitude - rightAscension;
		
		double c_lat = Math.cos(ConstPara.observerLatitude);
		double s_lat = Math.sin(ConstPara.observerLatitude);
		double c_H = Math.cos(hourAngle);
		double s_H = Math.sin(hourAngle);
		//parallax correction to right ascension:
		double d_alpha = -4.26e-5*c_lat*s_H;
		//topocentric right ascension
		topocRightAscension = rightAscension + d_alpha;
		//ropocentric hour angle
		topocHourAngle = hourAngle-d_alpha;
		//parallax correction to declination:
		topocDeclination = declination - 4.26e-5*(s_lat-declination*c_lat);
		double s_delta_corr = Math.sin(topocDeclination);
		double c_delta_corr = Math.cos(topocDeclination);
		double c_H_corr = c_H + d_alpha*s_H;
		double s_H_corr = s_H - d_alpha*c_H;
		
		//solar elevation angle,without refraction correction;
		elevationNoRefrac = Math.asin(s_lat*s_delta_corr + c_lat*c_delta_corr*c_H_corr);
		//refraction correction:it is calculated only if elevationNoRefrac>elev_min
		double elev_min = -0.01;
		if(elevationNoRefrac>elev_min) {
			//refractionCorrection = 0.084217*pressure/(273+temperature)/Math.tan(elevationNoRefrac + 0.0031376/(elevationNoRefrac+0.089186));
			refractionCorrection = 0.084217*ConstPara.pressure/
				((273+ConstPara.TA[period])*Math.tan(elevationNoRefrac + 0.0031376/ //wy131127
				 (elevationNoRefrac+0.089186)));
		} else {
			refractionCorrection = 0;
		}
		//local coordinates of the sun:
		zenith = Math.PI/2 - elevationNoRefrac - refractionCorrection;
		altitude = Math.PI/2 - zenith;
		azimuth = Math.atan2(s_H_corr, c_H_corr*s_lat-s_delta_corr/c_delta_corr*c_lat);	
		
		//energy of the sun， W/m^2
		int esc = 1367; //constant energy of the sun
		double esdcc; // earth-sun distance correction coefficient
		//esdcc = 1 + 0.0334*Math.cos(0.9856*whichDay(dd,mm,yy)-2.7206);
		esdcc = 1; //wy131127 ignore
		double sinh = Math.sin(altitude);
		float eos = (float)(esc * esdcc * sinh); //energy of the sun
		
		float directLight; //direct light
		float diffusedLight; //diffused light 
		if(eos<0) {
			directLight = 0.001f;
			diffusedLight = 0.1f;
		} else {
			directLight = (float)(eos*0.85); // 0.85 is the part of direct light in total radition 
			diffusedLight = eos - directLight;
		}
		
		//result 
		return float result[] = {
			(float)(altitude*180/Math.PI), 
			(float)(azimuth*180/Math.PI), 
			(float)(directLight), 
			(float)(diffusedLight)};
	}	

}