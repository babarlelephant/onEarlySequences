package beast.evolution.tree.coalescent;

public class modifiedCoalescent extends Coalescent {
	// fractional error less than x.xx * 10 ^ -4.
	// Algorithm 26.2.17 in Abromowitz and Stegun, Handbook of Mathematical.
	public static double erf2(double z) { // computes the primitive of exp(-x^2) vanishing at 0
		double t = 1.0 / (1.0 + 0.47047 * Math.abs(z));
		double poly = t * (0.3480242 + t * (-0.0958798 + t * (0.7478556)));
		double ans = 1.0 - poly * Math.exp(-z*z);
		if (z >= 0) return  ans;
		else        return -ans;
	}
	public static double Phi(double z) {
		return 1-erf2(Math.abs(z)); // computes  int_{|x|>|z|} exp(-x^2)dx
	}
	// number of cases by onset from December 8 to January 23
	int [] onsetCurve = new int[] {1,0,0,1,4,0,0,4,1,7,3,2,12,6,11,13,12,13,7,26,11,16,24,16,73,63,65,46,78,64,71,
					125,117,187,168,214,281,296,404,458,572,653,667,1254,1222,1552,1937};
							
	double serialInterval = 5/365.0; // 5 days
	int day0 = 32+25;  // January 25 (latest sequence in the dataset)
	double ascertainmentRate = 0.15; // our onset curve contains 0.15 less cases than the real number of infectees
	double sigma = 0.002;

	@Override
	public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction, double threshold) {
		ExponentialGrowth eg = (ExponentialGrowth) popSizeFunction;
		double logL = 0;
		for (int j = 0; j < onsetCurve.length; j++) {
			int day = j+8;
			int daysInPast = day-day0; // compared to day0 = Jan 25 
			double t = -daysInPast/365.0; // in year
			double effectivePopSize = eg.getPopSize(t);
			double realPopulationSize = effectivePopSize/serialInterval;
			double normalizedDiff = (onsetCurve[j]/ascertainmentRate -realPopulationSize) / Math.sqrt(realPopulationSize) * sigma;
			logL += Math.log(Phi(normalizedDiff));
		}
		return logL+super.calculateLogLikelihood(intervals,popSizeFunction,threshold);
	}
}
