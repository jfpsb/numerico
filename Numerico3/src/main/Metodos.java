package main;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class Metodos {

	private static List<Ponto> pontos = new ArrayList<Ponto>();

	public static void main(String[] args) {
		Newton.newton(new Array2DRowRealMatrix(new double[] { -1.2, 1 }), 0.0001, 0.0001);

		pontos.add(new Ponto(-1, 4));
		pontos.add(new Ponto(0, 1));
		pontos.add(new Ponto(2, -1));

		Lagrange.interpolar(17);

		InterporlacaoNewton.interpolar(17);

	}

	private static class Newton {
		public static void newton(RealMatrix x, double aproximacao1, double aproximacao2) {
			System.out.println("NEWTON:");

			RealMatrix F;
			RealMatrix jacobiano;
			RealMatrix s;

			while (true) {
				F = F(x);
				jacobiano = J(x);

				if (moduloMaximo(F) < aproximacao1) {
					printMatrix(x);
					break;
				} else {
					// Inversa do jacobiano
					RealMatrix inversa = new LUDecomposition(jacobiano).getSolver().getInverse();

					s = inversa.copy().multiply(F.copy().scalarMultiply(-1));

					x = x.add(s);

					if (moduloMaximo(s) < aproximacao2) {
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
			for (int i = 0; i < x.getRowDimension(); i++) {
				for (int j = 0; j < x.getColumnDimension(); j++) {
					System.out.print(x.getEntry(i, j));
					System.out.print("  ");
				}

				System.out.println();
			}
		}
	}

	private static class Lagrange {
		public static void interpolar(double x) {
			double y = 0;

			for (int i = 0; i < pontos.size(); i++) {
				y += (pontos.get(i).Y * L(x, i));
			}

			System.out.println("INTERPOLAÇÃO DE LAGRANGE: ");
			System.out.print(y);
		}

		public static double L(double x, int ordem) {
			double result = 1;

			for (int i = 0; i < pontos.size(); i++) {
				if (ordem != i) {
					result *= ((x - pontos.get(i).X) / (pontos.get(ordem).X - pontos.get(i).X));
				}
			}

			return result;
		}
	}

	private static class InterporlacaoNewton {
		public static void interpolar(double x) {
			double y = 0;
			double fator;

			for (int i = 0; i < pontos.size(); i++) {
				
				fator = 1;
				
				Ponto[] pontosDif = new Ponto[i + 1];
				
				for(int j = 0; j < i + 1; j++) {
					pontosDif[j] = pontos.get(j);
				}
				
				double diferenca = difDivididas(pontosDif);
				
				for(int w = 0; w < i; w++) {
					fator *= (x - pontos.get(w).X);
				}
				
				y += (fator * diferenca);
			}

			System.out.println("\nINTERPOLAÇÃO DE NEWTON: ");
			System.out.print(y);
		}
		
		private static double difDivididas(Ponto... pontos) {
			int ordem = pontos.length - 1;

			if (ordem == 0) {
				return pontos[0].Y;
			} else {
				Ponto[] pontos1 = new Ponto[pontos.length - 1];
				Ponto[] pontos2 = new Ponto[pontos.length - 1];

				for (int i = 1; i < pontos.length; i++) {
					pontos1[i - 1] = pontos[i];
				}

				for (int i = 0; i < pontos.length - 1; i++) {
					pontos2[i] = pontos[i];
				}

				return (difDivididas(pontos1) - difDivididas(pontos2)) / ((pontos[pontos.length - 1].X - pontos[0].X));
			}
		}
	}

	private static class Ponto {
		public double X;
		public double Y;

		public Ponto(double x, double y) {
			X = x;
			Y = y;
		}
	}
}
