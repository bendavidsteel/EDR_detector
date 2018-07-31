import java.util.Random;

public class Gene {

	double value;
	
	double lowerBound;
	double upperBound;
	
	Random rand = new Random();
	
	public Gene(double lowerBound, double upperBound) {
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		value = (upperBound - lowerBound)*rand.nextDouble() + lowerBound;
	}
	
	public void mutate() {
		value = (upperBound - lowerBound)*rand.nextDouble() + lowerBound;
	}
}
