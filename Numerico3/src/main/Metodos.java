package main;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fraction.Fraction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class Metodos {

	private static List<Ponto> pontos = new ArrayList<Ponto>();

	public static void main(String[] args) {
		
		EliminacaoDeGauss.chamarResolucao();
		
		// NÃO LINEAR NEWTON
		Newton.newton(new Array2DRowRealMatrix(new double[] { -1.2, 1 }), 0.0001, 0.0001);

		pontos.add(new Ponto(-1, 4));
		pontos.add(new Ponto(0, 1));
		pontos.add(new Ponto(2, -1));

		// INTERPOLAÇÃO LAGRANGE
		Lagrange.interpolar(17);

		//INTERPOLAÇÃO NEWTON
		InterporlacaoNewton.interpolar(17);
	}
	
	private static class EliminacaoDeGauss {
		
		public static void chamarResolucao() {
			System.out.println(">>> PÁGINA 373, QUESTÃO 8");
			System.out.println("LETRA A");
			System.out.println("SISTEMA:");
			
			RealMatrix letraAcoeficientes = new Array2DRowRealMatrix(3, 3);
			letraAcoeficientes.addToEntry(0, 0, 0.5);
			letraAcoeficientes.addToEntry(0, 1,0.25);
			letraAcoeficientes.addToEntry(0, 2, -0.125);
			letraAcoeficientes.addToEntry(1, 0, 0.33333333);
			letraAcoeficientes.addToEntry(1, 1, -0.16666666666);
			letraAcoeficientes.addToEntry(1, 2, 0.11111111111);
			letraAcoeficientes.addToEntry(2, 0, 0.14285714);
			letraAcoeficientes.addToEntry(2, 1, 0.14285714);
			letraAcoeficientes.addToEntry(2, 2, 0.1);
			
			RealMatrix letraAconstantes = new Array2DRowRealMatrix(3, 1);		
			letraAconstantes.addToEntry(0, 0, 0);
			letraAconstantes.addToEntry(1, 0, 1);
			letraAconstantes.addToEntry(2, 0, 2);
			
			printSistema(letraAcoeficientes, letraAconstantes);
			resolverSistema(letraAcoeficientes, letraAconstantes);
			
			System.out.println("LETRA B");
			System.out.println("SISTEMA:");
			
			RealMatrix letraBcoeficientes = new Array2DRowRealMatrix(3, 3);
			letraBcoeficientes.addToEntry(0, 0, 2.71);
			letraBcoeficientes.addToEntry(0, 1, 1);
			letraBcoeficientes.addToEntry(0, 2, 1032);
			letraBcoeficientes.addToEntry(1, 0, 4.12);
			letraBcoeficientes.addToEntry(1, 1, -1);
			letraBcoeficientes.addToEntry(1, 2, 500);
			letraBcoeficientes.addToEntry(2, 0, 3.33);
			letraBcoeficientes.addToEntry(2, 1, 2);
			letraBcoeficientes.addToEntry(2, 2, -200);
			
			RealMatrix letraBconstantes = new Array2DRowRealMatrix(3, 1);		
			letraBconstantes.addToEntry(0, 0, 12);
			letraBconstantes.addToEntry(1, 0, 11.49);
			letraBconstantes.addToEntry(2, 0, 41);
			
			printSistema(letraBcoeficientes, letraBconstantes);
			resolverSistema(letraBcoeficientes, letraBconstantes);
			
			System.out.println("LETRA C");
			System.out.println("SISTEMA:");
			
			RealMatrix letraCcoeficientes = new Array2DRowRealMatrix(4, 4);
			letraCcoeficientes.addToEntry(0, 0, 3.141592);
			letraCcoeficientes.addToEntry(0, 1, 1.414213);
			letraCcoeficientes.addToEntry(0, 2, -1);
			letraCcoeficientes.addToEntry(0, 3, 1);
			letraCcoeficientes.addToEntry(1, 0, 2.71828182846);
			letraCcoeficientes.addToEntry(1, 1, -1);
			letraCcoeficientes.addToEntry(1, 2, 1);
			letraCcoeficientes.addToEntry(1, 3, 2);
			letraCcoeficientes.addToEntry(2, 0, 1);
			letraCcoeficientes.addToEntry(2, 1, 1);
			letraCcoeficientes.addToEntry(2, 2, -1.73205080757);
			letraCcoeficientes.addToEntry(2, 3, 1);			
			letraCcoeficientes.addToEntry(3, 0, -1);
			letraCcoeficientes.addToEntry(3, 1, -1);
			letraCcoeficientes.addToEntry(3, 2, 1);
			letraCcoeficientes.addToEntry(3, 3, -2.2360679775);
			
			RealMatrix letraCconstantes = new Array2DRowRealMatrix(4, 1);		
			letraCconstantes.addToEntry(0, 0, 1);
			letraCconstantes.addToEntry(1, 0, 2);
			letraCconstantes.addToEntry(2, 0, 3);
			letraCconstantes.addToEntry(3, 0, 4);
			
			printSistema(letraCcoeficientes, letraCconstantes);
			resolverSistema(letraCcoeficientes, letraCconstantes);
			
			System.out.println("LETRA D");
			System.out.println("SISTEMA:");
			
			RealMatrix letraDcoeficientes = new Array2DRowRealMatrix(5, 5);
			letraDcoeficientes.addToEntry(0, 0, 1);
			letraDcoeficientes.addToEntry(0, 1, 1);
			letraDcoeficientes.addToEntry(0, 2, -1);
			letraDcoeficientes.addToEntry(0, 3, 1);
			letraDcoeficientes.addToEntry(0, 4, -1);
			letraDcoeficientes.addToEntry(1, 0, 2);
			letraDcoeficientes.addToEntry(1, 1, 2);
			letraDcoeficientes.addToEntry(1, 2, 1);
			letraDcoeficientes.addToEntry(1, 3, -1);
			letraDcoeficientes.addToEntry(1, 4, 1);
			letraDcoeficientes.addToEntry(2, 0, 3);
			letraDcoeficientes.addToEntry(2, 1, 1);
			letraDcoeficientes.addToEntry(2, 2, -3);
			letraDcoeficientes.addToEntry(2, 3, -2);
			letraDcoeficientes.addToEntry(2, 4, 3);
			letraDcoeficientes.addToEntry(3, 0, 4);
			letraDcoeficientes.addToEntry(3, 1, 1);
			letraDcoeficientes.addToEntry(3, 2, -1);
			letraDcoeficientes.addToEntry(3, 3, 4);
			letraDcoeficientes.addToEntry(3, 4, -5);
			letraDcoeficientes.addToEntry(4, 0, 16);
			letraDcoeficientes.addToEntry(4, 1, -1);
			letraDcoeficientes.addToEntry(4, 2, 1);
			letraDcoeficientes.addToEntry(4, 3, -1);
			letraDcoeficientes.addToEntry(4, 4, -1);
			
			RealMatrix letraDconstantes = new Array2DRowRealMatrix(5, 1);		
			letraDconstantes.addToEntry(0, 0, 2);
			letraDconstantes.addToEntry(1, 0, 4);
			letraDconstantes.addToEntry(2, 0, 8);
			letraDconstantes.addToEntry(3, 0, 16);
			letraDconstantes.addToEntry(4, 0, 32);
			
			printSistema(letraDcoeficientes, letraDconstantes);
			resolverSistema(letraDcoeficientes, letraDconstantes);
		}
		
		private static void printSistema(RealMatrix coeficientes, RealMatrix constantes) {
			for(int i = 0; i < coeficientes.getRowDimension(); i++) {
				for(int j = 0; j < coeficientes.getColumnDimension(); j++) {
					System.out.print(coeficientes.getEntry(i, j) + "x" + (j + 1) + "  ");
				}
				System.out.print("= " + constantes.getEntry(i, 0));
				System.out.println();
			}
		}
		
		public static void resolverSistema(RealMatrix coeficientes, RealMatrix constantes) {
			double pivo = 0;
			double multiplicador = 0;
			for(int i = 0; i < coeficientes.getRowDimension() - 1; i++) {
				pivo = coeficientes.getEntry(i, i);
				
				if(pivo == 0) {
					for(int p = i + 1; p < coeficientes.getRowDimension(); p++) {
						double pivo_candidato = coeficientes.getEntry(p, i);
						
						if(pivo_candidato != 0) {
							double[] linhaAtual;
							double[] linhaPivoCandidato;
							
							linhaAtual = coeficientes.getRow(i).clone();
							linhaPivoCandidato = coeficientes.getRow(p).clone();
							
							coeficientes.setRow(i, linhaPivoCandidato);
							coeficientes.setRow(p, linhaAtual);
							
							pivo = coeficientes.getEntry(i, i);
							
							break;
						}
					}
				}
				
				for(int j = i + 1; j < coeficientes.getRowDimension(); j++) {
					multiplicador = coeficientes.getEntry(j, i)/pivo;
					
					for(int w = 0; w < coeficientes.getColumnDimension(); w++) {
						double resultCoeficientes = coeficientes.getEntry(j, w) - (multiplicador * coeficientes.getEntry(i, w));
						coeficientes.setEntry(j, w, resultCoeficientes);
					}
					
					double resultConstantes = constantes.getEntry(j, 0) - (multiplicador * constantes.getEntry(i, 0));
					constantes.setEntry(j, 0, resultConstantes);
				}
			}
			
			System.out.println("\nMATRIZ TRIANGULAR SUPERIOR:");
			for(int i = 0; i < coeficientes.getRowDimension(); i++) {
				for (int j = 0; j < coeficientes.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", coeficientes.getEntry(i, j));
				}
				System.out.printf("| %.6f", constantes.getEntry(i, 0));
				System.out.println();
			}
			
			resolveTriagularSuperior(coeficientes, constantes);
		}
		
		public static void resolveTriagularSuperior(RealMatrix coeficientes, RealMatrix constantes) {
			double soma = 0;
			double variavel = 0;
			
			RealMatrix variaveis = new Array2DRowRealMatrix(coeficientes.getRowDimension(), 1);
			
			variaveis.setEntry(coeficientes.getRowDimension() - 1, 0, constantes.getEntry(coeficientes.getRowDimension() - 1, 0)/coeficientes.getEntry(coeficientes.getRowDimension() - 1, coeficientes.getRowDimension() - 1));
			
			for(int i = coeficientes.getRowDimension() - 2; i >= 0; i--) {
				soma = 0;
				
				// Retrosubstituição
				for(int j = i + 1; j < coeficientes.getRowDimension(); j++) {
					soma += coeficientes.getEntry(i, j) * variaveis.getEntry(j, 0);
				}
				
				variavel = (constantes.getEntry(i, 0) - soma)/(coeficientes.getEntry(i, i));
				
				variaveis.addToEntry(i, 0, variavel);
			}
			
			System.out.println("\nSOLUCÃO:");
			
			for(int i = 0; i < variaveis.getRowDimension(); i++) {
				System.out.printf("x%d = %.6f\n", i, variaveis.getEntry(i, 0));
			}
			
			System.out.println();
		}
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
