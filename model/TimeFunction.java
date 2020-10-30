/**
 * @descr update time
 * @author weiyang 20131127
 */
public class TimeFunction {
	
	//update date of year
	public static int[] updateTimeByDay(int year, int month, int day) {
		int[] dom = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};//days of month
		if( (year%400==0) || (year%4 == 0 && year%100!=0) ) {
			dom[1] = 29;
		}
		//change day of month
		if(month == 12) {
			if(day>=1 && day<dom[month-1]) {                                                                                                            
				day ++;
			} else if(day == dom[month-1]) {
				day = 1;
				month = 1;
				year ++;
			}
		} else {
			if(day>=1 && day<dom[month-1]) {
				day ++;
			} else if(day == dom[month-1]) {
				day = 1;
				month ++;
			}
		}
		return int result[] = {year, month, day};
	}
	
	//calculate the date is which day of the year
	public static int whichDayOfTheYear(int year, int month, int day)
	{
		int wd = 0;
		int[] dayBase = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334}; //day base of each month
		//leap year
		if( (year%400==0) || (year%4 == 0 && year%100!=0) )
		{
			dayBase[2] += 1;
		}
		wd = dayBase[month-1] + day;
		return wd;
	}

}