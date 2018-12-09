package main;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class Metodos {

	public static void main(String[] args) {
		Newton(new Array2DRowRealMatrix(new double[] { -1.2, 1 }), 0.0001, 0.0001);

	}

	public static void Newton(RealMatrix x, double aproximacao1, double aproximacao2) {
		RealMatrix F;
		RealMatrix jacobiano;
		RealMatrix s;

		while (true) {
			F = F(x);
			jacobiano = J(x);

			if (moduloMaximo(F) < aproximacao1) {
				System.out.print("RESULTADO: ");
				printMatrix(x);
				break;
			} else {
				// Inversa do jacobiano
				RealMatrix inversa = new LUDecomposition(jacobiano).getSolver().getInverse();

				s = inversa.copy().multiply(F.copy().scalarMultiply(-1));

				x = x.add(s);

				if (moduloMaximo(s) < aproximacao2) {
					System.out.print("RESULTADO: ");
					printMatrix(x);
					break;
				}
			}
		}
	}

	public static RealMatrix F(RealMatrix x) {
		Array2DRowRealMatrix result = new Array2DRowRealMatrix(2, 1);

		result.setEntry(0, 0, eq1(x.getEntry(0, 0), x.getEntry(1, 0)));
		result.setEntry(1, 0, eq2(x.getEntry(0, 0), x.getEntry(1, 0)));

		return result;
	}

	public static double eq1(double x1, double x2) {
		return (10 * x2) - (10 * Math.pow(x1, 2));
	}

	public static double eq2(double x1, double x2) {
		return 1 - x1;
	}

	public static RealMatrix J(RealMatrix x) {
		RealMatrix j = new Array2DRowRealMatrix(2, 2);

		j.setEntry(0, 0, (20 * x.getEntry(0, 0)) * (-1));
		j.setEntry(0, 1, 10);
		j.setEntry(1, 0, -1);
		j.setEntry(1, 1, 0);

		return j;
	}

	public static double moduloMaximo(RealMatrix F) {
		double maior = Math.abs(F.getEntry(0, 0));

		if (Math.abs(F.getEntry(1, 0)) > maior) {
			maior = Math.abs(F.getEntry(1, 0));
		}

		return maior;
	}

	public static void printMatrix(RealMatrix x) {
		System.out.println();
		
		for (int i = 0; i < x.getRowDimension(); i++) {
			for (int j = 0; j < x.getColumnDimension(); j++) {
				System.out.print(x.getEntry(i, j));
				System.out.print("  ");
			}
			
			System.out.println();
		}
	}

}
