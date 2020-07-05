/*

Java-adaptation of the following LEAFC3-model:
 
(C) 2001 by Ned Nikolov, Ph.D.

    N & T Services, LLC
    6736 Fragler Rd.
    Fort Collins, CO 80525
    Phone: (980) 980-3303
    E-mail: ntconsulting@comcast.net

Revised version of the Delphi Pascal code of the LEAFC3 model published in:


   Nikolov, N.T., Massman,W.J. & Schoettle, A.W. 1995. Coupling biochemical
   and biophysical processes at the leaf level: an equilibrium photosynthesis
   model for leaves of C3 plants. Ecol. Modelling 80:205-235.


The interface section of the unit assigns default values to a number
of model parameters (see declaration const in the implementation section) that
do not normally change with every call of procedure LeafC3. These parameters may
vary as a function of plant species and/or season and environmental conditions.
Thus, they have been declared in a way which allows programs using this unit to
assign new values to them whenever appropriate. The parameters CCgs and CCgb are
convergence criteria for Gs and Gb, respectively. They have been
assigned default values that optimize both computational speed and
numerical accuracy.

For more details about input parameters and output variables, and instructions
of how to use this unit see Nikolov et al. 1995.
*/


/*
Environmental input parameters:

Mm     Soil moisture multiplier (0 <= Mm <=1.0) It quantifies the effect of soil 
       drought on stomatal conductance.
Pa     Atmospheric pressure,                  Pa
Ca     Ambient CO2 concentration,             mol/mol
Oa     Ambient oxygen concentration,          mol/mol
Ta     Ambient air temperature,               deg C
RH     Ambient relative humidity,             decimal fraction
Tdew   Dew-point air temperature              deg C
Ri     Bi-directional absorbed total short-and long-wave radiation by the leaf, 
       W/m2
PPFD   Incident photosynthetic photon flux density, umol m-2 s-1
WSpeed Horizontal Wind speed,                 m/s
Wstat  Leaf wetness status,                   0 - dry; 1 - wet
ShFac  Shelter factor accounting for the effect of foliage aggregation on leaf 
       boundary-layer conductance inside the canopy (0 < ShFac <= 1). Typically, 
       ShFac = 0.5;

NB. The parameter Ri is NOT to be confused with the absorbed net radiation!
    In this model:
    Ri = (1-Wv)*Qv + (1-Wi)*Qi + (LWd + LWu)
where Wv and Wi are leaf scattering coefficients (reflectance + transmittance)
in the visible and infrared band, respectively; Qv and Qi are incident visible
and infrared radiative fluxes (W m-2); LWd and LWu are long-wave (thermal)
radiative fluxes (W m-2) coming from above and below, respectively.
    PPFD = 4.6*Qv          or
    PPFD = 2.162*(Qv + Qi)
The latter equation is an approximation assuming that global radiation flux
consists of 50% visible and 50% near-infrared flux.

 */
 
class Environment
{
	public Environment (float IR) {this.IR=IR;}
	float IR;
	//final Parameter par = Parameter.INSTANCE; //wy131125
	final float latitude = ConstPara.observerLatitude2;
	public double Mm = 1.0; 
	public double Alt = 10; //GBS100430 /* altitude above sea level, value for Los Banos, Philippines*/
	//public double Pa = 101325.0; 
	public double Ca = 0.000375; /* Ambient CO2 concentration, mol/mol */
	public double Oa = 0.209; /* Ambient oxygen concentration, mol/mol */
	public double[] Ta = ConstPara.TA;
	/* Ambient air temperature, °C, 20 cm above soil.
	 * Values from IRRI wetland weather station,  Los Baños, Laguna
	 * , 1.3. - 30.6.2008 */
	//GBS100501
	public double SimTempCourse (double Temp, int h)
	{
		double F = 6.0;
		double S = 16.4;
		double Freq = 3.6;
		return (Temp + F * Math.sin((h/Freq) + S));
		//return Temp;
	} 
	/* simulated daily course of temperature; sine function based on average temperature */

	//public double[] dawn = par.DAWN;
	//public double[] dusk = par.DUSK; 
	public double SimPPFDCourse (int doy, int time) //GBS100501
	{
		double corrTime;
		double ecc = 0.01671; //eccentricity
		double ax = 23.438; //angle of Earth axis
		float daymn1 = 2.0f;
		float daymn2 = 81.0f;
		double rad = Math.PI/180;
		double corr1 = -4.0*2.0*ecc*Math.sin(2*Math.PI*(doy-daymn1)/365.0f)/rad;
		double corr2 = (ax*rad)*(ax*rad)*Math.sin(4*Math.PI*(doy-daymn2)/365.0f)/rad;
		double corr = corr1 + corr2;
		corrTime = time + corr;
		double sinl = Math.sin(rad * latitude); //ok
		double cosl = Math.cos(rad * latitude); //ok
		//sine of declination, sine delta:
		double sind = -Math.sin(ax*rad)*Math.cos(2*Math.PI*(doy+10)/365); //ok
		//cosine of declination, cosine delta:
		double cosd = Math.sqrt(1-(sind*sind)); //ok
		//cosine of hour angle:
		//float costh = Math.cos(rad*15*(time-12)); //ok
		double costh = Math.cos(2*Math.PI*(time-12)/24);
		//sine of solar height, sine alpha
		double sina = (sinl * sind) + (cosl * cosd * costh);
		double S = Math.max(0,grtemp[doy]* sina); //GBS100322
		//double S = Math.max(0,1367* sina); //GBS100322
		return S; 
	}
	/* simulated daily course of PPFD; sinus function based on maximum global radiation at Midday */
	
	//public double[] tempsum = par.TEMP_SUM; //GBS100501: not used here
	/* Temperature sum, °C d, 20 cm above soil. 
	 * Values from IRRI wetland weather station,  Los Baños, Laguna,
	 * 1.3. - 30.6.2008 */
	
	//public double[] daylength = par.DAY_LENGTH; //GBS100501: not used here
	
	public double Pa (int t) 
	{
		return 101325.0*Math.exp(-Alt/(29.3*(Ta[t] + 273.16)));
	}
	/* Atmospheric pressure, Pa, standard value taken, modified by altitude and temperature */
	public double[] RH = ConstPara.REL_HUMIDITY;  
	/* Ambient relative humidity. 
	 * Values from IRRI wetland weather station,  Los Baños, Laguna,
	 * 1.3. - 30.6.2008 */
	public double Tdew (int t, int h)
	{
	 //return  237.3/((7.5/((Math.log(RH[t]*100)/2.30258509)+(7.5*Ta[t]/(237.3+Ta[t]))-2))-1); 
	 return  273.3/((7.5/((Math.log(RH[t]/* GBS100501 */*100)/2.30258509)+(7.5*SimTempCourse(Ta[t],h)/(273.3+SimTempCourse(Ta[t],h)))-2))-1);
	}
	/* Dew point temperature  */
	public double[] gr = {IR};//par.gr;
	/* Global radiation (Qv + Qi). 
	 * Values from IRRI wetland weather station,  Los Baños, Laguna,
	 * 1.3. - 30.6.2008 */
	
	public double[] grtemp = {556, 563, 571, 579, 586, 594, 602, 610, 617, 625, 633, 640, 648, 
			656, 663, 671, 678, 686, 693, 701, 708, 715, 723, 730, 737, 744, 751, 758, 765, 
			772, 779, 786, 793, 799, 806, 813, 819, 825, 832, 838, 844, 850, 856, 862, 868, 
			874, 880, 885, 891, 896, 901, 907, 912, 917, 922, 927, 932, 937, 941, 946, 950, 
			955, 959, 963, 967, 971, 975, 979, 983, 987, 990, 994, 997, 1001, 1004, 1007, 
			1010, 1013, 1016, 1019, 1021, 1024, 1027, 1029, 1031, 1034, 1036, 1038, 1040, 
			1042, 1044, 1045, 1047, 1049, 1050, 1052, 1053, 1054, 1055, 1057, 1058, 1059, 
			1059, 1060, 1061, 1062, 1062, 1063, 1063, 1063, 1063, 1064, 1064, 1064, 1064, 
			1063, 1063, 1063, 1062, 1062, 1061, 1061, 556, 563, 571, 579, 586, 594, 602, 610, 617, 625, 633, 640, 648, 
			656, 663, 671, 678, 686, 693, 701, 708, 715, 723, 730, 737, 744, 751, 758, 765, 
			772, 779, 786, 793, 799, 806, 813, 819, 825, 832, 838, 844, 850, 856, 862, 868, 
			874, 880, 885, 891, 896, 901, 907, 912, 917, 922, 927, 932, 937, 941, 946, 950, 
			955, 959, 963, 967, 971, 975, 979, 983, 987, 990, 994, 997, 1001, 1004, 1007, 
			1010, 1013, 1016, 1019, 1021, 1024, 1027, 1029, 1031, 1034, 1036, 1038, 1040, 
			1042, 1044, 1045, 1047, 1049, 1050, 1052, 1053, 1054, 1055, 1057, 1058, 1059, 
			1059, 1060, 1061, 1062, 1062, 1063, 1063, 1063, 1063, 1064, 1064, 1064, 1064, 
			1063, 1063, 1063, 1062, 1062, 1061, 1061}; 
	public double Qv (int t)
	{
		//return gr[t]*0.5;
		return IR;
	}/*  incoming solar radiation in the visible band (W m-2) */
	public double PPFD (int t)
	{
		return Qv (t) * 4.6; 
	}
	/* Incident photosynthetic photon flux density, µmol m-2 s-1 */
	//public double Qi = PPFD*0.5 / 2.162; /*  incoming radiation in the near-infrared band (W m-2) */
	//public double Cloud = 0.5; /* cloud cover */
	public double[] Cloud = ConstPara.CLOUD; 
	public double Esoil = 0.965; /* thermal emissivity of the ground (the soil) */
	public double Ea = 1200; /* atmospheric water-vapor pressure (Pa) */
	public double atmosLWRad (double Ta, double Ea, double Cloud, int type)
	/* Returns incident long-wave atmospheric radiation (W m-2) using screen-height 
	 * weather data.
	 * Input:
	 * Ta = air temperature (°C)
	 * Ea = atmos. water-vapor pressure (Pa)
	 * Cloud = fractional cloud cover (0 <= Cloud <= 1) */
	{
	double emis, AtmosLWRad;
	if (type == 1)
	{
		emis = (1 - 0.38*Math.exp(-10.6e-2*Ea/(Ta+273.16)))*(1 + 0.005*Math.pow(8*Cloud,2));
		AtmosLWRad = 5.67e-8*emis*Math.pow(Math.pow(Ta + 273.16,2),2);
	}
	else 
	AtmosLWRad = 5.67e-8*Esoil*Math.pow(Math.pow(Ta + 273.16,2),2);
	return AtmosLWRad;
	}
	
	public double Wv = 0.13;
	public double Wi = 0.88;
	public double Ri (int t) 
	{
		return (1-Wv)*Qv(t) + (1-Wi)*Qv(t) + (atmosLWRad(Ta[t],Ea,Cloud[t],1) + atmosLWRad(Ta[t],Ea,Cloud[t],0));
	}
	/* Bi-directional absorbed short- and long-wave radiation by the leaf, W/m2. */
	public double[] WSpeed = ConstPara.WSPEED;  
	/*  Wind velocity, m/s. 
	 *  Values from IRRI wetland weather station,  Los Baños, Laguna,
	 *  1.3. - 30.6.2008 */
	public double Wstat = 0; /* Leaf wetness status, 0 - dry, 1 - wet */
	public double ShFac = 0.9 /* 0.5 */;
}

/*
Species-specific input parameters:

Vm25   Maximum carboxylation velocity at 25C,  umol m-2 s-1
Jm25   Light-saturated potential rate of electron transport at 25C, µmol m-2 s-1
Ej     Activation energy for electron transport J/mol
Hj     Enthalpy parameter,                     J/mol
Sj     Enthalpy parameter,                     J mol-1 K-1
theta  Coefficient controling the smoothness of transition between light and
    temperature limitations on the potential rate of electron transport --
Cdr    Proportionality coefficient for estimating leaf night respiration rate
    from Vm25;
Kc25   Kinetic parameter for CO2 at 25C,       mol/mol
Ko25   Kinetic parameter for O2 at 25C,        mol/mol
f      PPFD loss factor,                       -
m      Composite stomatal sensitivity,         -
Bs     Empirical constant for stomatal conductance, mol m-2 s-1
Dleaf  Leaf width (or needle diameter),        m
Dshoot Shoot diameter (for conifers only),     m
*/
class Species
{
	public double Vm25 = 118.9; /* µmol m-2 s-1, value for rice, Yin et al., 2004, p. 1217, Tab. 5 "New model" */
	public double Jm25 = 194.5; /* µmol m-2 s-1, value for rice, Yin et al., 2004, p. 1217, Tab. 5 "New model" */
	
	public double Ej = 888380; /* Activation energy, J/mol, value for rice, Yin 2005, Table 2 and 3, pp. 46-47 */
	public double Hj = 219.7e3; /* 219.7e3; *//* Enthalpy parameter, J/mol */
	public double Sj /* = 710.0; */ = 650; /* Enthalpy parameter, J mol-1 K-1, from Harley et al. 1992 */
	
	public double Kc25 = 0.00040409; /* mol/mol, value for wheat, Müller et al. 2005, Table C.1, p. 208 */
	public double Ko25 = 0.2784; /* mol/mol, value for wheat, Müller et al. 2005, Table C.1, p. 208  */
	public double f = 0.4; /* -, Müller, p. 204 */
	public double m = 7.34; /* -, value for mature leaves, Müller, p. 204 */
	public double Bs = 0.0376; /* mol m-2 s-1, Müller, Eq. (9) */
	public double Dleaf = 0.01; /* m,  value for rice, Yin 2005, Table 2 and 3, pp. 46-47 */
	public double Dshoot = 0.01; /* m, estimated value for cereals */

	public double Aql = 0.87; /* Leaf absorption coeficient for PPFD */

	public boolean Astom = false; /* True = amphystomatous leaf, False = hypostomatous leaf */

    /* Coefficient controling the smoothness
            of transition between light and
            temperature limitations on the potential
            rate of electron transport */
	public double theta = 0.8;

	/* Proportionality coefficient for
	 	estimating leaf night respiration rate */
	public double Cdr = 0.0089; /* Watanabe, Evans, & Chow 1994, value for T = 25°C */
	

	/* Returns leaf net photosynthesis rate (umol m-2 s-1) by solving
	   a quartic equation. The function also returns leaf maintenance (dark)
	   respiration (Rd) (µmol m-2 s-1).*/
	double C3NetAssim
		(final double[] AgrOut, final double[] RdOut, final double Pa,
		 final double Cap, final double Oa, final double PPFD,
		 final double Iss, final double Gs, final double Gb,
		 final double Tl, final double Mbv)
	{
		double Km,Ccp,Vcm,Jcm,J,R,a,sa,b,c,d,L,alfa,gm,Rt,sRt,N,Q,P,W,t1,VmQ10;
		double Wc,Wj;
		
		N = Tl + 273.16;
	    R = 0.12027905/N;
	    L = (PPFD < 50) ? PPFD*0.01332 : 0.666;
	    double Rd = Cdr*Vm25*(1.666 - L)*Math.exp(0.074193734*(Tl-25)); // assumes Rd Q10 = 2.1
	    double Agr;
	    if (Tl > 49)
	    {
	    	Rd /= (1 + Math.exp(1.3*(Tl-55)));
	    }
	    if ((PPFD > 0) && (Gs > 1e-3))
	    {
	    	Q = Tl - 25;
	    	Ccp = Pa*Oa*(213.88e-6 + Q*(8.995e-6 + 1.722e-7*Q));
	    	Km = Pa*Kc25*Math.exp(32.462 - 80470*R);
	    	if (Tl <= 0)
	    	{
	    		VmQ10 = 0.13083328;	// Q10 = 3.7
	    	}
	    	else if (Tl < 25.0)
	    	{
	    		VmQ10 = 0.1*Math.log(3.7 - 0.055*Tl); // variable Q10
	    	}
	        else
	        {
	        	VmQ10 = 0.08501509; // Q10 = 2.34
	        }
	    	Vcm = Vm25*Math.exp(VmQ10*(Tl-25));
	    	if (Tl > 28)
	    	{
	    		Vcm /= (1 + Math.exp(0.9*(Tl-41.5)));
	    	}
	    	Jcm = Jm25*Math.exp((N*3.362135e-3 - 1)*Ej*R)/(1 + Math.exp((Sj*N - Hj)*R));
	    	J = 0.111111*(Jcm + Iss - Math.sqrt((Jcm+Iss)*(Jcm+Iss) - 4*theta*Jcm*Iss))/theta;
	    	Rt = Pa*(1.577*Mbv*Gb + 1.355*Gs)/(Mbv*Gs*Gb);
	    	sRt = Rt*Rt;
	    	Km = Km*(1 + Oa/(Ko25*Math.exp(5.854 - 14510/(8.314*N))));
	    	N = (2.33*Ccp + Km);
	    	P = Vcm + J;
	    	W = Vcm*J;
	    	Q = 1.33*Vcm*Ccp + J*(Km - Ccp);
	    	L = Ccp*(2.33*Vcm*Ccp + J*Km);
	    	alfa = 0.96*(Cap*(Cap+N) + 2.33*Ccp*Km);
	    	gm = 1/(0.96*sRt);
	    	R = 0.96*(2*Cap+N);
	    	a = gm*(sRt*(2*0.96*Rd - P) - Rt*R);
	    	b = gm*(alfa + Rt*(Rd*(Rt*(0.96*Rd - P) - 2*R) +
	    		2*P*Cap + Q + W*Rt));
	    	c = gm*(Rd*(2*alfa + Rt*(2*P*Cap + Q - Rd*R)) -
	    		Cap*(P*Cap + Q) + L - 2*W*Rt*(Cap - Ccp));
	    	d = gm*(Rd*(Rd*alfa - Cap*(P*Cap + Q) + L) + W*(Cap*(Cap - 2*Ccp) +
	    		Ccp*Ccp));
	    	sa = a*a;
	    	gm = a*c - 4*d;
	    	L = b*b;
	    	P = Math.abs((3*gm - L))*0.1111111111;
	    	W = Math.sqrt(P);
	    	Q = (b*(2*L-9*gm) - 27*(d*(4*b-sa)-c*c))/54;
	    	t1 = Q/(P*W);
	    	t1 = Math.acos(t1)*0.3333333333;
	    	N = 2*W*Math.cos(t1) + 0.3333333333*b;
	    	R = Math.sqrt(0.25*sa + N - b);
	    	L = Math.sqrt(0.5*sa - b - N - 0.25*(a*(4*b-sa)-8*c)/R);
	    	Agr = -0.25*a - 0.5*(R+L) + Rd;
	    	Wc = 0.37*Vcm;
	    	Agr = 0.515464*(Agr+Wc - Math.sqrt((Agr + Wc)*(Agr + Wc) - 3.88*Agr*Wc));
	    }
	    else
	    {
	    	Agr = 0;
	    }
	    AgrOut[0] = Agr;
	    RdOut[0] = Rd;
	    return Agr - Rd;
	}

}


public class LeafC3
{
	/*
	Output Variables

	An     Net photosynthesis rate,              umol m-2 s-1
	Resp   Leaf maintenance respiration,         umol m-2 s-1
	Gs     Stomatal conductance to water vapor,  mmol m-2 s-1
	Gb     All-sided leaf-boundary layer conductance to water vapor, mmol m-2 s-1
	GtO3   Total leaf conductance to ozone exchange,       mm/s
	GtCO2  Total leaf conductance to CO2 exchange,       mmol m-2 s-1
	Tl     Leaf temperature,                     deg C
	LHeat  Latent heat flux from the leaf,       W/m^2
	*/

	public double An;
	public double Resp;
	public double Gs;
	public double Gb;
	public double GtO3;
	public double GtCO2;
	public double Tl;
	public double LHeat;
	public double Td;
	public int d = 122; /* max. Day in weather array: 1..121 */
	//public int i;

    /* Convergence criterium for stomatal
            conductance, umol m-2 s-1 */
	public static final double CCgs = 100;
	
    /* Convergence criterium for leaf-boundary
            layer conductance, umol m-2 s-1 */
	public static final double CCgb = 3000;

	/* wy131120add
	 * This is the main program of the LEAFC3 model
	 * @para env: environment
	 * @para spec: species
	 * @para i: day
	 * @para h: hour
	 */
	public void compute (final Environment env, final Species spec, int i, int h)
	{
		double Agr,Mt,LHV,Cfm,G,Ea,Es,Bsc,rPs,rP,dGb,Gbf,Grf,Gbm,Gm,
			VHCair,Gbc,Mbv,Tav,Iss,Cap,Ria,Tlt,Est;
		int n;
		double[] sc = new double[6];
		boolean FreeConvection;
		
		Cap = env.Pa(i)*env.Ca; 
		//SVPCoef (sc,env.Ta[i],'w');
		SVPCoef (sc,env.SimTempCourse(env.Ta[i],h),'w'); //GBS100501
		//Es = (env.Ta[i]*(sc[4]+env.Ta[i]*(sc[3]+env.Ta[i]*(sc[2]+sc[1]*
		//	env.Ta[i]))) + sc[5]);
		//GBS100501:
		Es = (env.SimTempCourse(env.Ta[i],h)*(sc[4]+env.SimTempCourse(env.Ta[i],h)*(sc[3]+env.SimTempCourse(env.Ta[i],h)*(sc[2]+sc[1]*
			env.SimTempCourse(env.Ta[i],h)))) + sc[5]);
		Ea = env.RH[i]*Es;
		rP = 1/env.Pa(i);
		//Tav = env.Ta[i] + 273.16;
		Tav = env.SimTempCourse(env.Ta[i],h) + 273.16; //GBS100501
		Cfm = 8.3089764e-6*Tav*rP;
		//LHV = (2.50084 - 0.00234*env.Ta[i])*1e6; // latent heat of vaporization, J/kg
		LHV = (2.50084 - 0.00234*env.SimTempCourse(env.Ta[i],h))*1e6; //GBS100501
		rPs = 2.1669236e-3*LHV/Tav;
		VHCair = 3.518638*env.Pa(i)/Tav;
		Bsc = (env.SimPPFDCourse(i,h) < 100) ? 1e6*spec.Bs*(env.SimPPFDCourse(i,h)*0.009 + 0.1) : 1e6*spec.Bs;
		Iss = 0.5*spec.Aql*(1-spec.f)*env.SimPPFDCourse(i,h);
		Mt = spec.m*env.Mm;
		Gbf = 1.6693*Math.sqrt(env.WSpeed[i]*env.Pa(i)/spec.Dleaf);
		if ((spec.Dshoot > 0) && (spec.Dleaf <= 5e-3))
		{
			Gbf = env.ShFac*144.84335*Gbf;
			Mbv = 1;
		}
		else
		{
			Gbf = env.ShFac*520.16034*Gbf;
			Mbv = spec.Astom ? 1 : 0.5;
		}
		Gb = Math.max (Gbf, 0.2e6);
		FreeConvection = (Gbf <= 0.2e6) || (env.WSpeed[i] <= 0.28);
		if (FreeConvection)
		{
			Tav = Tav/(1 - 0.378*Ea*rP);
			Grf = ((spec.Dshoot > 0) && (spec.Dleaf <= 5e-3))
				? 174.162*Math.sqrt(env.Pa(i)/Math.sqrt(spec.Dshoot))
				: 328.697*Math.sqrt(env.Pa(i)/Math.sqrt(spec.Dleaf));
		}
		else
		{
			Grf = 0;
		}
		Gbm = Gb*Cfm;
		
		// Estimates initial value for Gs
		/* Gs = (env.SimPPFDCourse(i,h) > 0)
			? (env.SimPPFDCourse(i,h)/(env.SimPPFDCourse(i,h)+65))*0.15*spec.Vm25*Math.exp(0.07*(env.Ta[i]-25))*Mt*env.RH[i]/env.Ca + Bsc
			: 0.5*Bsc; */
		//GBS100501:
		Gs = (env.SimPPFDCourse(i,h) > 0)
			? (env.SimPPFDCourse(i,h)/(env.SimPPFDCourse(i,h)+65))*0.15*spec.Vm25*Math.exp(0.07*(env.SimTempCourse(env.Ta[i],h)-25))*Mt*env.RH[i]/env.Ca + Bsc
			: 0.5*Bsc; 
		n = 0;
		dGb = 0;
		Ria = env.Ri(i);
		//Tl = env.Ta[i];
		Tl = env.SimTempCourse(env.Ta[i],h); //GBS100501
		final double[] var0 = new double[1], var1 = new double[1];
		do
		{
			n++;
			G = Gs;
			Gm = G*Cfm;
			if (FreeConvection)
			{
				Gbc = LBLCond (Gbf,Grf,Tl,env.SimTempCourse(env.Ta[i],h),Tav,Ea,Es,rP,Gs,env.Wstat); //GBS100501
				//Gbc = LBLCond (Gbf,Grf,Tl,env.Ta[i],Tav,Ea,Es,rP,Gs,env.Wstat);
				dGb = Gb - Gbc;
				Gb = Gbc;
				Gbm = Gb*Cfm;
			}
			Tl = SurfTemp (var0,Ria,env.SimTempCourse(env.Ta[i],h),Ea,rPs,VHCair,Gm,Gbm,sc,-1,0,0,env.Wstat,Mbv,0.975); //GBS100501
			//Tl = SurfTemp (var0,Ria,env.Ta[i],Ea,rPs,VHCair,Gm,Gbm,sc,-1,0,0,env.Wstat,Mbv,0.975);
			Es = var0[0];
			if ((Tl - env.Tdew(i,h)) < 2.2)
			//if ((Tl - env.Tdew(i)) < 2.2)
			{
				Tlt = SurfTemp (var0,Ria,env.SimTempCourse(env.Ta[i],h),Ea,rPs,VHCair,Gm,Gbm,sc,-1,0,0,1,Mbv,0.975); //GBS100501
				Est = var0[0];
				if (Tlt < env.Tdew(i,h))
				//if (Tlt < env.Tdew(i))
				{
					Tl = Tlt;
					Es = Est;
				}
			}
			An = spec.C3NetAssim (var0,var1,env.Pa(i),Cap,env.Oa,/* env.SimPPFDCourse(i,h)*/env.PPFD(i),Iss,G,Gb,Tl,Mbv);
			//System.out.println("LEAFC3 - An: " + An);
			Agr = var0[0];
			Resp = var1[0];
			Gs = (An > 0) ? StomCond (Es,Mt,Bsc,An,env.Ca,Gb,Ea,env.Wstat,Mbv) : Bsc;
			Ria = env.Ri(i) - 0.506*Agr; /* Accounts for the energy consumed in
					photosynthetic reactions, i.e. 0.506 J/umol CO2.*/
		} while (!(((Math.abs(Gs-G) < CCgs) && (Math.abs(dGb) < CCgb)) || (n > 45)));
		LHeat = rPs*(Es - Ea)*Cfm*Gb;
		if ((env.Wstat <= 0) && (LHeat > 0))
		{
			LHeat *= Gs*Mbv/(Gs + Mbv*Gb);
		}
		G = Gs*Gb*Mbv;
		GtO3 = 1e3*Cfm*G/(1.315*Gs + Mbv*1.508*Gb);
		
		    // GtO3 = total leaf conductance to ozone exchange in mm/s
		GtCO2 = 1e-3*G/(1.355*Gs + Mbv*1.577*Gb);
		    // GtCO2 = total leaf conductance to CO2 exchange in mmol m-2 s-1
		Td = env.Tdew(i,h);
		//Td = env.Tdew(i); //GBS100501
		Gs *= 1e-3;
		Gb *= 1e-3;
		}
	
	/* Assignes 5 coefficients (array sc) of a fourth-order polynomial used
			   to estimate the saturation vapor pressure (SVP) over water (base = 'w')
			   or ice (base = 'i'). The general form of the polynomial is:
			     SVP = sc[1]*T^4 + sc[2]*T^3 + sc[3]*T^2 + sc[4]*T + sc[5]
			   where T is air temperature (C), and SVP is estimated in Pascal.
			   The returned coefficients allow for an accurate computation of SVP
			   over the range -50 C to +50 C. The polynomial is used in the analytical
			   solution of the energy balance equation for surface temperature (see function
			   SurfTemp below).
			Input:
			   Temp = ambient air temperature (deg C)
			   base = 'w' (water) or 'i' (ice).
			*/
	public static void SVPCoef
		(double[] sc, final double Temp, final char base)
	{
		if (base != 'i')
		{
			if (Temp > 9)
			{
				sc[1] = 6.0810032e-4f;
				sc[2] = 1.3979668e-2f;
				sc[3] = 1.5979711f;
				sc[4] = 44.27433f;
				sc[5] = 608.56329f;
			}
		    else if (Temp > -7.0)
		    {
		    	sc[1] = 3.0242495e-4f;
		    	sc[2] = 2.7794843e-2f;
		    	sc[3] = 1.4323892f;
		    	sc[4] = 44.281468f;
		    	sc[5] = 611.38265f;
		    }
		    else
		    {
		    	sc[1] = 1.1560863e-4f;
		    	sc[2] = 1.9467751e-2f;
		    	sc[3] = 1.3193403f;
		    	sc[4] = 43.836964f;
		    	sc[5] = 611.21804f;
		    }
		}
		else
		{
			sc[1] = 1.4921906e-4f;
			sc[2] = 2.4687197e-2f;
			sc[3] = 1.5820521f;
			sc[4] = 48.041683f;
			sc[5] = 601.733f;
		}
	}
	
	
	/*
			 Returns surface temperature (deg C) and saturation vapor pressure (Es, in Pa)
			 at the surface by analytically solving a quartic form of the energy-balance
			 equation (no loops involved!). The function is applicable to foliage, soil and
			 snow surfaces.
			INPUT:
			 Ri = total absorbed radiation (W m-2) by the surface. For foliage, Ri refers
			      to the radiation absorbed from two directions. For soil and snow, Ri is
			      the absorbed part of radiation incident from above only.
			 Ta = air temperature outside the surface boundary layer (C)
			 Ea = ambient atmos. water vapor pressure (Pascals)
			 rPs = Cp*Da/Cps, where Cp is specific heat of dry air (J kg-1 K-1), Da is
			       density of air (kg m-3), Cps is the psychrometric constant (Pa/K).
			 HCair = volumetric heat capacity of air (J m-3 K-1)
			 Gs = surface conductance to water vapor (m/s)
			 Gb = boundary layer conductance to water vapor (m/s)
			 sc = an array containing coefficents of a fourth order polynomial used to
			      estimating saturation vapor pressure (see procedure SVPCoef above).
			 Wstat = surface wetness status (0-dry; 1-wet)
			 Mbv = parameter specifying hypostomatal (0.5) or amphystomatal (1.0)
			       leaf morphology (applicable to foliage surfaces only)
			 Kt = thermal conductivity of the medium (W m-1 K-1). Note, for leaves, Kt
			      should have a value equal to or less than zero.
			 Tz = temperature of the medium at depth Z (deg C). This parameter is only
			      used when estimating soil- or snowpack-surface temperature.
			 Z =  depth (m) at which Tz is measured.
			 emis = thermal emissivity of the surface.
			*/
	public static double SurfTemp
		(final double[] Es, final double Ri, final double Ta,
		 final double Ea, final double rPs, final double HCair,
		 final double Gs, final double Gb,
		 double[] sc, final double Kt, double Tz, final double Z,
		 final double Wstat, final double Mbv, final double emis)
	{
		double he,ht,k,k1,a,sa,b,c,d,Q,P,Dscr,y,R,t1,E;
		double Ts;

		if (Kt > 0)
		{
			k1 = Kt/Z;
			t1 = 0.5*emis;
		}
		else
		{
			k1 = 0;
			Tz = 0;
			t1 = emis;
		}
		he = (Wstat <= 0) ? rPs*Gs*Gb*Mbv/(Gs + Mbv*Gb) : rPs*Gb;
		ht = HCair*(0.924*Gb);
	    k = 1/(t1*1.10565e-7 + he*sc[1]);
	    a = k*(t1*1.2080774e-4 + he*sc[2]);
	    b = k*(t1*0.049499764 + he*sc[3]);
	    c = k*(t1*9.014237 + he*sc[4] + ht + k1);
	    d = k*(he*(sc[5]-Ea) - Ri - k1*Tz - ht*Ta + t1*615.58224);
	    y = a*c - 4*d;
	    E = b*b;
	    sa = a*a;
	    P = (3*y - E)*0.1111111111;
	    Q = (b*(2*E - 9*y) - 27*(d*(4*b-sa)-c*c))/54;
	    Dscr = Math.sqrt(Q*Q + P*P*P);
	    y = Math.exp(Math.log(Q+Dscr)/3) - Math.exp(Math.log(Dscr-Q)/3) + 0.3333333333*b;
	    R = Math.sqrt(0.25*sa + y - b);
	    E = 0.5*sa - b - y - 0.25*(a*(4*b-sa)-8*c)/R;
	    Ts = -0.25*a - 0.5*(R-Math.sqrt(E));
	    Es[0] = sc[5] + Ts*(sc[4] + Ts*(sc[3] + Ts*(sc[2] + sc[1]*Ts)));
	    return Ts;
	}


	/*
			  Returns leaf-boundary layer conductance (umol m-2 s-1) as the larger of
			  the conductances resulting from forced- and free-convective exchanges.
			*/

	private static double LBLCond
		(final double Gbf, final double Grf, final double Tl,
		 final double Ta, final double Tav, final double Ea,
		 final double Es, final double rP, final double Gs,
		 final double Wstat)
	{
		double Gb, Eb, Gbe, Tlk;
		int n;

		Tlk = Tl + 273.16;
	    if (Wstat <= 0)
	    {
	    	Eb = Math.abs(Tl-Ta);
	    	Gb = (Eb > 1) ? Grf*Math.sqrt(Math.sqrt(Eb*(1 + 0.378*rP*Ea))) : 0.2e6;
	    	n = 0;
	    	do
	    	{
	    		n++;
	    		Gbe = Gb;
	    		Eb = (Gs*Es + Gbe*Ea)/(Gbe+Gs);
	    		Eb = Tlk/(1 - 0.378*Eb*rP) - Tav;
	    		Gb = Grf*Math.sqrt(Math.sqrt(Math.abs(Eb)));
	    	} while ((Math.abs(Gb-Gbe) > CCgb) && (n < 12));
	    }
	    else
	    {
	    	Eb = Tlk/(1 - 0.378*Es*rP) - Tav;
	    	Gb = Grf*Math.sqrt(Math.sqrt(Math.abs(Eb)));
	    }
	    return Math.max (Gbf, Gb);
	}


	/*
			 Calculates leaf stomatal conductance to water vapor (umol m-2 s-1) by
			 solving a quadratic equation.
			*/
	private static double StomCond
		(final double Es, final double M, final double Bs,
		 final double An, final double Ca, final double Gb,
		 final double Ea, final double Wstat, final double Mbv)
	{
		double Cb, a, b, c, MAn, Gs, Gbt;
		
		MAn = M*An;
		Gbt = Mbv*Gb;
		Cb = Ca - 1.355*An/Gbt;
		if (Wstat <= 0)
		{
			a = Cb;
			b = MAn - Cb*(Gbt-Bs);
			c = Gbt*(MAn*Ea/Es + Cb*Bs);
			Gs = 0.5*(b + Math.sqrt(b*b + 4*a*c))/a;
		}
		else
		{
			Gs = MAn/Cb + Bs;
		}
		return Math.max (Gs, 0);
	}
	
}
