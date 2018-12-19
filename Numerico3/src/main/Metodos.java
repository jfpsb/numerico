package main;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class Metodos {

	private static List<Ponto> pontos = new ArrayList<Ponto>();

	public static void main(String[] args) {

		 GaussSeidel.chamarResolucao();

		Choleski.chamarResolucao();

		 FatoracaoLU.chamarResolucao();

		 EliminacaoDeGauss.chamarResolucao();

		// NÃO LINEAR NEWTON
		// Newton.newton(new Array2DRowRealMatrix(new double[] { -1.2, 1 }), 0.0001,
		// 0.0001);

		// pontos.add(new Ponto(-1, 4));
		// pontos.add(new Ponto(0, 1));
		// pontos.add(new Ponto(2, -1));

		// INTERPOLAÇÃO LAGRANGE
		// Lagrange.interpolar(17);

		// INTERPOLAÇÃO NEWTON
		// InterporlacaoNewton.interpolar(17);
	}

	private static class GaussJordan {
		public static void chamarResolucao() {
			RealMatrix coeficientes = new Array2DRowRealMatrix(3, 3);
			RealMatrix constantes = new Array2DRowRealMatrix(3, 1);

			coeficientes.setEntry(0, 0, 3);
			coeficientes.setEntry(0, 1, -1);
			coeficientes.setEntry(0, 2, 1);
			coeficientes.setEntry(1, 0, 3);
			coeficientes.setEntry(1, 1, 6);
			coeficientes.setEntry(1, 2, 2);
			coeficientes.setEntry(2, 0, 3);
			coeficientes.setEntry(2, 1, 3);
			coeficientes.setEntry(2, 2, 7);

			constantes.setEntry(0, 0, 1);
			constantes.setEntry(1, 0, 0);
			constantes.setEntry(2, 0, 4);

			GaussJordan.resolverSistema(coeficientes, constantes);
		}
		
		public static void resolverSistema(RealMatrix coeficientes, RealMatrix constantes) {
			double pivo = 0;
			RealMatrix coeficientes2 = coeficientes.copy();
			RealMatrix constantes2 = constantes.copy();
			
			for (int i = 0; i < coeficientes.getRowDimension(); i++) {
				pivo = coeficientes.getEntry(i, i);
				for (int j = 0; j < coeficientes.getColumnDimension(); j++) {
					coeficientes.setEntry(i, j, coeficientes.getEntry(i, j) / pivo);
				}

				constantes.setEntry(i, 0, constantes.getEntry(i, 0) / pivo);
				
				coeficientes2.setRow(i, coeficientes.getRow(i));
				constantes2.setColumn(0, constantes.getColumn(0));

				for (int p = 0; p < coeficientes.getRowDimension(); p++) {
					if (p != i) {
						for (int y = 0; y < coeficientes.getColumnDimension(); y++) {
							double aux = coeficientes2.getEntry(p, y) - (coeficientes2.getEntry(p, i) * coeficientes2.getEntry(i, y));
							coeficientes.setEntry(p, y, aux);
							constantes.setEntry(p, 0, constantes2.getEntry(p, 0) - (coeficientes2.getEntry(p, i) * constantes2.getEntry(i, 0)));
						}
					}
				}
			}
		}
	}

	private static class GaussSeidel {
		public static void chamarResolucao() {
			System.out.println("GAUSS SEIDEL");
			System.out.println("PÁGINA 465, QUESTÃO 1");
			System.out.println("LETRA A");
			RealMatrix letraACoeficientes = new Array2DRowRealMatrix(3, 3);
			RealMatrix letraAConstantes = new Array2DRowRealMatrix(3, 1);
			RealMatrix aproximacao = new Array2DRowRealMatrix(3, 1);
			RealMatrix aproximacao2 = new Array2DRowRealMatrix(4, 1);
			RealMatrix aproximacao3 = new Array2DRowRealMatrix(5, 1);

			aproximacao.setEntry(0, 0, 0);
			aproximacao.setEntry(1, 0, 0);
			aproximacao.setEntry(2, 0, 0);

			aproximacao2.setEntry(0, 0, 0);
			aproximacao2.setEntry(1, 0, 0);
			aproximacao2.setEntry(2, 0, 0);
			aproximacao2.setEntry(3, 0, 0);

			aproximacao3.setEntry(0, 0, 0);
			aproximacao3.setEntry(1, 0, 0);
			aproximacao3.setEntry(2, 0, 0);
			aproximacao3.setEntry(3, 0, 0);
			aproximacao3.setEntry(4, 0, 0);

			letraACoeficientes.setEntry(0, 0, 3);
			letraACoeficientes.setEntry(0, 1, -1);
			letraACoeficientes.setEntry(0, 2, 1);

			letraACoeficientes.setEntry(1, 0, 3);
			letraACoeficientes.setEntry(1, 1, 6);
			letraACoeficientes.setEntry(1, 2, 2);

			letraACoeficientes.setEntry(2, 0, 3);
			letraACoeficientes.setEntry(2, 1, 3);
			letraACoeficientes.setEntry(2, 2, 7);

			letraAConstantes.setEntry(0, 0, 1);
			letraAConstantes.setEntry(1, 0, 0);
			letraAConstantes.setEntry(2, 0, 4);

			GaussSeidel.resolverSistema(letraACoeficientes, aproximacao, letraAConstantes, 0.00001);

			System.out.println("\nLETRA B");
			RealMatrix letraBCoeficientes = new Array2DRowRealMatrix(3, 3);
			RealMatrix letraBConstantes = new Array2DRowRealMatrix(3, 1);

			letraBCoeficientes.setEntry(0, 0, 10);
			letraBCoeficientes.setEntry(0, 1, -1);
			letraBCoeficientes.setEntry(0, 2, 0);

			letraBCoeficientes.setEntry(1, 0, -1);
			letraBCoeficientes.setEntry(1, 1, 10);
			letraBCoeficientes.setEntry(1, 2, -2);

			letraBCoeficientes.setEntry(2, 0, 0);
			letraBCoeficientes.setEntry(2, 1, -2);
			letraBCoeficientes.setEntry(2, 2, 10);

			letraBConstantes.setEntry(0, 0, 9);
			letraBConstantes.setEntry(1, 0, 7);
			letraBConstantes.setEntry(2, 0, 6);

			GaussSeidel.resolverSistema(letraBCoeficientes, aproximacao, letraBConstantes, 0.00001);

			System.out.println("\nLETRA C");
			RealMatrix letraCCoeficientes = new Array2DRowRealMatrix(4, 4);
			RealMatrix letraCConstantes = new Array2DRowRealMatrix(4, 1);

			letraCCoeficientes.setEntry(0, 0, 10);
			letraCCoeficientes.setEntry(0, 1, 5);
			letraCCoeficientes.setEntry(0, 2, 0);
			letraCCoeficientes.setEntry(0, 3, 0);

			letraCCoeficientes.setEntry(1, 0, 5);
			letraCCoeficientes.setEntry(1, 1, 10);
			letraCCoeficientes.setEntry(1, 2, -4);
			letraCCoeficientes.setEntry(1, 3, 0);

			letraCCoeficientes.setEntry(2, 0, 0);
			letraCCoeficientes.setEntry(2, 1, -4);
			letraCCoeficientes.setEntry(2, 2, 8);
			letraCCoeficientes.setEntry(2, 3, -1);

			letraCCoeficientes.setEntry(3, 0, 0);
			letraCCoeficientes.setEntry(3, 1, 0);
			letraCCoeficientes.setEntry(3, 2, -1);
			letraCCoeficientes.setEntry(3, 3, 5);

			letraCConstantes.setEntry(0, 0, 6);
			letraCConstantes.setEntry(1, 0, 25);
			letraCConstantes.setEntry(2, 0, -11);
			letraCConstantes.setEntry(3, 0, -11);

			GaussSeidel.resolverSistema(letraCCoeficientes, aproximacao2, letraCConstantes, 0.00001);

			System.out.println("\nLETRA D");
			RealMatrix letraDCoeficientes = new Array2DRowRealMatrix(5, 5);
			RealMatrix letraDConstantes = new Array2DRowRealMatrix(5, 1);

			letraDCoeficientes.setEntry(0, 0, 4);
			letraDCoeficientes.setEntry(0, 1, 1);
			letraDCoeficientes.setEntry(0, 2, 1);
			letraDCoeficientes.setEntry(0, 3, 0);
			letraDCoeficientes.setEntry(0, 4, 1);

			letraDCoeficientes.setEntry(1, 0, -1);
			letraDCoeficientes.setEntry(1, 1, -3);
			letraDCoeficientes.setEntry(1, 2, 1);
			letraDCoeficientes.setEntry(1, 3, 1);
			letraDCoeficientes.setEntry(1, 4, 0);

			letraDCoeficientes.setEntry(2, 0, 2);
			letraDCoeficientes.setEntry(2, 1, 1);
			letraDCoeficientes.setEntry(2, 2, 5);
			letraDCoeficientes.setEntry(2, 3, -1);
			letraDCoeficientes.setEntry(2, 4, -1);

			letraDCoeficientes.setEntry(3, 0, -1);
			letraDCoeficientes.setEntry(3, 1, -1);
			letraDCoeficientes.setEntry(3, 2, -1);
			letraDCoeficientes.setEntry(3, 3, 4);
			letraDCoeficientes.setEntry(3, 4, 0);

			letraDCoeficientes.setEntry(4, 0, 0);
			letraDCoeficientes.setEntry(4, 1, 2);
			letraDCoeficientes.setEntry(4, 2, -1);
			letraDCoeficientes.setEntry(4, 3, 1);
			letraDCoeficientes.setEntry(4, 4, 4);

			letraDConstantes.setEntry(0, 0, 6);
			letraDConstantes.setEntry(1, 0, 6);
			letraDConstantes.setEntry(2, 0, 6);
			letraDConstantes.setEntry(3, 0, 6);
			letraDConstantes.setEntry(4, 0, 6);

			GaussSeidel.resolverSistema(letraDCoeficientes, aproximacao3, letraDConstantes, 0.00001);
		}

		public static void resolverSistema(RealMatrix coeficientes, RealMatrix aproximacao, RealMatrix constantes,
				double precisao) {
			double soma = 0;
			RealMatrix aproximacaoAnterior;
			RealMatrix aproximacaoDiferenca = new Array2DRowRealMatrix(aproximacao.getRowDimension(), 1);

			while (true) {
				aproximacaoAnterior = aproximacao.copy();

				for (int i = 0; i < aproximacao.getRowDimension(); i++) {
					soma = constantes.getEntry(i, 0);
					for (int j = 0; j < coeficientes.getColumnDimension(); j++) {
						if (i != j)
							soma -= coeficientes.getEntry(i, j) * aproximacao.getEntry(j, 0);
					}

					soma /= coeficientes.getEntry(i, i);

					aproximacao.setEntry(i, 0, soma);
				}

				for (int i = 0; i < aproximacao.getRowDimension(); i++) {
					aproximacaoDiferenca.setEntry(i, 0,
							Math.abs(aproximacaoAnterior.getEntry(i, 0) - aproximacao.getEntry(i, 0)));
				}

				double d = Newton.moduloMaximo(aproximacaoDiferenca) / Newton.moduloMaximo(aproximacao);

				if (d < precisao) {
					System.out.println("SOLUÇÃO:");
					Newton.printMatrix(aproximacao);
					break;
				}
			}
		}
	}

	private static class Choleski {
		public static void chamarResolucao() {
			System.out.println("CHOLESKI");
			System.out.println("PÁGINA 429, QUESTÃO 3");
			System.out.println("LETRA A");
			RealMatrix letraA = new Array2DRowRealMatrix(3, 3);

			letraA.setEntry(0, 0, 2);
			letraA.setEntry(0, 1, -1);
			letraA.setEntry(0, 2, 0);

			letraA.setEntry(1, 0, -1);
			letraA.setEntry(1, 1, 2);
			letraA.setEntry(1, 2, -1);

			letraA.setEntry(2, 0, 0);
			letraA.setEntry(2, 1, -1);
			letraA.setEntry(2, 2, 2);

			fatorar(letraA);
			fatoracaoLDL(letraA);

			System.out.println("\nLETRA B");
			RealMatrix letraB = new Array2DRowRealMatrix(4, 4);

			letraB.setEntry(0, 0, 4);
			letraB.setEntry(0, 1, 1);
			letraB.setEntry(0, 2, 1);
			letraB.setEntry(0, 3, 1);

			letraB.setEntry(1, 0, 1);
			letraB.setEntry(1, 1, 3);
			letraB.setEntry(1, 2, -1);
			letraB.setEntry(1, 3, 1);

			letraB.setEntry(2, 0, 1);
			letraB.setEntry(2, 1, -1);
			letraB.setEntry(2, 2, 2);
			letraB.setEntry(2, 3, 0);

			letraB.setEntry(3, 0, 1);
			letraB.setEntry(3, 1, 1);
			letraB.setEntry(3, 2, 0);
			letraB.setEntry(3, 3, 2);

			fatorar(letraB);
			fatoracaoLDL(letraB);

			System.out.println("\nLETRA C");
			RealMatrix letraC = new Array2DRowRealMatrix(4, 4);

			letraC.setEntry(0, 0, 4);
			letraC.setEntry(0, 1, 1);
			letraC.setEntry(0, 2, -1);
			letraC.setEntry(0, 3, 0);

			letraC.setEntry(1, 0, 1);
			letraC.setEntry(1, 1, 3);
			letraC.setEntry(1, 2, -1);
			letraC.setEntry(1, 3, 0);

			letraC.setEntry(2, 0, 1);
			letraC.setEntry(2, 1, -1);
			letraC.setEntry(2, 2, 5);
			letraC.setEntry(2, 3, 2);

			letraC.setEntry(3, 0, 0);
			letraC.setEntry(3, 1, 0);
			letraC.setEntry(3, 2, 2);
			letraC.setEntry(3, 3, 4);

			fatorar(letraC);
			fatoracaoLDL(letraC);

			System.out.println("\nLETRA D");
			RealMatrix letraD = new Array2DRowRealMatrix(4, 4);

			letraD.setEntry(0, 0, 6);
			letraD.setEntry(0, 1, 2);
			letraD.setEntry(0, 2, 1);
			letraD.setEntry(0, 3, -1);

			letraD.setEntry(1, 0, 2);
			letraD.setEntry(1, 1, 4);
			letraD.setEntry(1, 2, 1);
			letraD.setEntry(1, 3, 0);

			letraD.setEntry(2, 0, 1);
			letraD.setEntry(2, 1, 1);
			letraD.setEntry(2, 2, 4);
			letraD.setEntry(2, 3, -1);

			letraD.setEntry(3, 0, -1);
			letraD.setEntry(3, 1, 0);
			letraD.setEntry(3, 2, -1);
			letraD.setEntry(3, 3, 3);

			fatorar(letraD);
			fatoracaoLDL(letraD);
		}

		public static void fatorar(RealMatrix coeficientes) {
			double soma;
			double soma2;
			RealMatrix G = new Array2DRowRealMatrix(coeficientes.getRowDimension(), coeficientes.getColumnDimension());

			System.out.println("\nFATORAÇÃO CHOLESKI");

			for (int k = 0; k < coeficientes.getRowDimension(); k++) {
				soma = 0;
				// Elementos da diagonal
				for (int j = 0; j <= k; j++) {
					soma += Math.pow(G.getEntry(k, j), 2);
				}

				G.setEntry(k, k, Math.sqrt(coeficientes.getEntry(k, k) - soma));

				// Outros elementos
				for (int i = k + 1; i < coeficientes.getRowDimension(); i++) {
					soma2 = 0;

					for (int j = 0; j <= k; j++) {
						soma2 += G.getEntry(i, j) * G.getEntry(k, j);
					}

					G.setEntry(i, k, (coeficientes.getEntry(i, k) - soma2) / G.getEntry(k, k));
				}
			}

			System.out.println("MATRIZ G:");
			for (int i = 0; i < G.getRowDimension(); i++) {
				for (int j = 0; j < G.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", G.getEntry(i, j));
				}
				System.out.println();
			}

			System.out.println("\nMATRIZ G^t:");
			RealMatrix Gt = G.transpose();
			for (int i = 0; i < Gt.getRowDimension(); i++) {
				for (int j = 0; j < Gt.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", Gt.getEntry(i, j));
				}
				System.out.println();
			}
		}

		public static void fatoracaoLDL(RealMatrix coeficientes) {
			double pivo = 0;
			double multiplicador = 0;

			System.out.println("\nFATORAÇÃO LDL^t");

			for (int i = 0; i < coeficientes.getRowDimension() - 1; i++) {
				pivo = coeficientes.getEntry(i, i);

				for (int j = i + 1; j < coeficientes.getRowDimension(); j++) {
					multiplicador = coeficientes.getEntry(j, i) / pivo;

					for (int w = i; w < coeficientes.getColumnDimension(); w++) {
						double resultCoeficientes = coeficientes.getEntry(j, w)
								- (multiplicador * coeficientes.getEntry(i, w));
						coeficientes.setEntry(j, w, resultCoeficientes);
					}

					coeficientes.setEntry(j, i, multiplicador);
				}
			}

			RealMatrix L = new Array2DRowRealMatrix(coeficientes.getRowDimension(), coeficientes.getColumnDimension());
			RealMatrix U_linha = new Array2DRowRealMatrix(coeficientes.getRowDimension(),
					coeficientes.getColumnDimension());
			RealMatrix D = new Array2DRowRealMatrix(coeficientes.getRowDimension(), coeficientes.getColumnDimension());

			// Popula matriz L
			for (int i = 0; i < L.getRowDimension(); i++) {
				L.setEntry(i, i, 1);

				for (int j = 0; j < L.getColumnDimension(); j++) {
					if (i > j) {
						L.setEntry(i, j, coeficientes.getEntry(i, j));
					}
				}
			}

			// Popula matriz D
			for (int i = 0; i < coeficientes.getRowDimension(); i++) {
				D.setEntry(i, i, coeficientes.getEntry(i, i));
			}

			// U é transposta de L
			U_linha = L.transpose();

			System.out.println("MATRIZ L:");
			for (int i = 0; i < L.getRowDimension(); i++) {
				for (int j = 0; j < L.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", L.getEntry(i, j));
				}
				System.out.println();
			}

			System.out.println("\nMATRIZ D:");
			for (int i = 0; i < D.getRowDimension(); i++) {
				for (int j = 0; j < D.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", D.getEntry(i, j));
				}
				System.out.println();
			}

			System.out.println("\nMATRIZ(L^t):");
			for (int i = 0; i < U_linha.getRowDimension(); i++) {
				for (int j = 0; j < U_linha.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", U_linha.getEntry(i, j));
				}
				System.out.println();
			}
		}
	}

	private static class FatoracaoLU {
		public static void chamarResolucao() {
			System.out.println("FATORAÇÃO LU");
			System.out.println(">>> PÁGINA 414, QUESTÃO 5");
			System.out.println("LETRA A");
			System.out.println("SISTEMA:");
			RealMatrix letraAcoeficientes = new Array2DRowRealMatrix(3, 3);
			letraAcoeficientes.addToEntry(0, 0, 2);
			letraAcoeficientes.addToEntry(0, 1, -1);
			letraAcoeficientes.addToEntry(0, 2, 1);
			letraAcoeficientes.addToEntry(1, 0, 3);
			letraAcoeficientes.addToEntry(1, 1, 3);
			letraAcoeficientes.addToEntry(1, 2, 9);
			letraAcoeficientes.addToEntry(2, 0, 3);
			letraAcoeficientes.addToEntry(2, 1, 3);
			letraAcoeficientes.addToEntry(2, 2, 5);
			
			resolverSistema(letraAcoeficientes, null);
			
			System.out.println("LETRA B");
			System.out.println("SISTEMA:");			
			RealMatrix letraBcoeficientes = new Array2DRowRealMatrix(3, 3);
			letraBcoeficientes.addToEntry(0, 0, 1.012);
			letraBcoeficientes.addToEntry(0, 1, -2.132);
			letraBcoeficientes.addToEntry(0, 2, 3.104);
			letraBcoeficientes.addToEntry(1, 0, -2.132);
			letraBcoeficientes.addToEntry(1, 1, 4.096);
			letraBcoeficientes.addToEntry(1, 2, -7.013);
			letraBcoeficientes.addToEntry(2, 0, 3.104);
			letraBcoeficientes.addToEntry(2, 1, -7.013);
			letraBcoeficientes.addToEntry(2, 2, 0.014);
			
			resolverSistema(letraBcoeficientes, null);
			
			System.out.println("LETRA C");
			System.out.println("SISTEMA:");			
			RealMatrix letraCcoeficientes = new Array2DRowRealMatrix(4, 4);
			letraCcoeficientes.addToEntry(0, 0, 2);
			letraCcoeficientes.addToEntry(0, 1, 0);
			letraCcoeficientes.addToEntry(0, 2, 0);
			letraCcoeficientes.addToEntry(0, 3, 0);
			letraCcoeficientes.addToEntry(1, 0, 1);
			letraCcoeficientes.addToEntry(1, 1, 1.5);
			letraCcoeficientes.addToEntry(1, 2, 0);
			letraCcoeficientes.addToEntry(1, 3, 0);
			letraCcoeficientes.addToEntry(2, 0, 0);
			letraCcoeficientes.addToEntry(2, 1, -3);
			letraCcoeficientes.addToEntry(2, 2, 0.5);
			letraCcoeficientes.addToEntry(2, 3, 0);
			letraCcoeficientes.addToEntry(3, 0, 2);
			letraCcoeficientes.addToEntry(3, 1, -2);
			letraCcoeficientes.addToEntry(3, 2, 1);
			letraCcoeficientes.addToEntry(3, 3, 1);
			
			resolverSistema(letraCcoeficientes, null);
			
			System.out.println("LETRA D");
			System.out.println("SISTEMA:");
			RealMatrix letraDcoeficientes = new Array2DRowRealMatrix(4, 4);
			letraDcoeficientes.addToEntry(0, 0, 2.1756);
			letraDcoeficientes.addToEntry(0, 1, 4.0231);
			letraDcoeficientes.addToEntry(0, 2, -2.1732);
			letraDcoeficientes.addToEntry(0, 3, 5.1967);
			letraDcoeficientes.addToEntry(1, 0, -4.0231);
			letraDcoeficientes.addToEntry(1, 1, 6);
			letraDcoeficientes.addToEntry(1, 2, 0);
			letraDcoeficientes.addToEntry(1, 3, 1.1973);
			letraDcoeficientes.addToEntry(2, 0, -1);
			letraDcoeficientes.addToEntry(2, 1, -5.2107);
			letraDcoeficientes.addToEntry(2, 2, 1.1111);
			letraDcoeficientes.addToEntry(2, 3, 0);
			letraDcoeficientes.addToEntry(3, 0, 6.0235);
			letraDcoeficientes.addToEntry(3, 1, 7);
			letraDcoeficientes.addToEntry(3, 2, 0);
			letraDcoeficientes.addToEntry(3, 3, -4.1561);
			
			resolverSistema(letraDcoeficientes, null);
		}

		public static void resolverSistema(RealMatrix coeficientes, RealMatrix constantes) {
			double pivo = 0;
			double multiplicador = 0;
			for (int i = 0; i < coeficientes.getRowDimension() - 1; i++) {
				pivo = coeficientes.getEntry(i, i);

				for (int j = i + 1; j < coeficientes.getRowDimension(); j++) {
					multiplicador = coeficientes.getEntry(j, i) / pivo;

					for (int w = i; w < coeficientes.getColumnDimension(); w++) {
						double resultCoeficientes = coeficientes.getEntry(j, w)
								- (multiplicador * coeficientes.getEntry(i, w));
						coeficientes.setEntry(j, w, resultCoeficientes);
					}

					coeficientes.setEntry(j, i, multiplicador);
				}
			}

			RealMatrix L = new Array2DRowRealMatrix(coeficientes.getRowDimension(), coeficientes.getColumnDimension());
			RealMatrix U = new Array2DRowRealMatrix(coeficientes.getRowDimension(), coeficientes.getColumnDimension());

			for (int i = 0; i < L.getRowDimension(); i++) {
				L.setEntry(i, i, 1);

				for (int j = 0; j < L.getColumnDimension(); j++) {
					if (i > j) {
						L.setEntry(i, j, coeficientes.getEntry(i, j));
					}
				}
			}

			for (int i = 0; i < U.getRowDimension(); i++) {
				for (int j = 0; j < U.getColumnDimension(); j++) {
					if (i <= j) {
						U.setEntry(i, j, coeficientes.getEntry(i, j));
					}
				}
			}

			System.out.println("\nMATRIZ L:");
			for (int i = 0; i < L.getRowDimension(); i++) {
				for (int j = 0; j < L.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", L.getEntry(i, j));
				}
				System.out.println();
			}

			System.out.println("\nMATRIZ U:");
			for (int i = 0; i < U.getRowDimension(); i++) {
				for (int j = 0; j < U.getColumnDimension(); j++) {
					System.out.printf("%.6f  ", U.getEntry(i, j));
				}
				System.out.println();
			}

//			System.out.println("RESOLVENDO Ly = b");
//			System.out.println("y:");
//			RealMatrix y = new LUDecomposition(L).getSolver().getInverse().multiply(constantes);
//			Newton.printMatrix(y);
//			
//			System.out.println("RESOLVENDO Ux = y");
//			System.out.println("x:");
//			RealMatrix x = new LUDecomposition(U).getSolver().getInverse().multiply(y);
//			Newton.printMatrix(x);
		}
	}

	private static class EliminacaoDeGauss {

		public static void chamarResolucao() {
			System.out.println("ELIMINAÇÃO DE GAUSS");
			System.out.println(">>> PÁGINA 373, QUESTÃO 8");
			System.out.println("LETRA A");
			System.out.println("SISTEMA:");

			RealMatrix letraAcoeficientes = new Array2DRowRealMatrix(3, 3);
			letraAcoeficientes.addToEntry(0, 0, 0.5);
			letraAcoeficientes.addToEntry(0, 1, 0.25);
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
			for (int i = 0; i < coeficientes.getRowDimension(); i++) {
				for (int j = 0; j < coeficientes.getColumnDimension(); j++) {
					System.out.print(coeficientes.getEntry(i, j) + "x" + (j + 1) + "  ");
				}
				System.out.print("= " + constantes.getEntry(i, 0));
				System.out.println();
			}
		}

		public static void resolverSistema(RealMatrix coeficientes, RealMatrix constantes) {
			double pivo = 0;
			double multiplicador = 0;
			for (int i = 0; i < coeficientes.getRowDimension() - 1; i++) {
				pivo = coeficientes.getEntry(i, i);

				if (pivo == 0) {
					for (int p = i + 1; p < coeficientes.getRowDimension(); p++) {
						double pivo_candidato = coeficientes.getEntry(p, i);

						if (pivo_candidato != 0) {
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

				for (int j = i + 1; j < coeficientes.getRowDimension(); j++) {
					multiplicador = coeficientes.getEntry(j, i) / pivo;

					for (int w = 0; w < coeficientes.getColumnDimension(); w++) {
						double resultCoeficientes = coeficientes.getEntry(j, w)
								- (multiplicador * coeficientes.getEntry(i, w));
						coeficientes.setEntry(j, w, resultCoeficientes);
					}

					double resultConstantes = constantes.getEntry(j, 0) - (multiplicador * constantes.getEntry(i, 0));
					constantes.setEntry(j, 0, resultConstantes);
				}
			}

			System.out.println("\nMATRIZ TRIANGULAR SUPERIOR:");
			for (int i = 0; i < coeficientes.getRowDimension(); i++) {
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

			variaveis.setEntry(coeficientes.getRowDimension() - 1, 0,
					constantes.getEntry(coeficientes.getRowDimension() - 1, 0) / coeficientes
							.getEntry(coeficientes.getRowDimension() - 1, coeficientes.getRowDimension() - 1));

			for (int i = coeficientes.getRowDimension() - 2; i >= 0; i--) {
				soma = 0;

				// Retrosubstituição
				for (int j = i + 1; j < coeficientes.getRowDimension(); j++) {
					soma += coeficientes.getEntry(i, j) * variaveis.getEntry(j, 0);
				}

				variavel = (constantes.getEntry(i, 0) - soma) / (coeficientes.getEntry(i, i));

				variaveis.addToEntry(i, 0, variavel);
			}

			System.out.println("\nSOLUCÃO:");

			for (int i = 0; i < variaveis.getRowDimension(); i++) {
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

			for (int i = 1; i < F.getRowDimension(); i++) {
				if (Math.abs(F.getEntry(i, 0)) > maior) {
					maior = Math.abs(F.getEntry(i, 0));
				}
			}

			return maior;
		}

		public static void printMatrix(RealMatrix x) {
			for (int i = 0; i < x.getRowDimension(); i++) {
				for (int j = 0; j < x.getColumnDimension(); j++) {
					System.out.printf("%.6f", x.getEntry(i, j));
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

				for (int j = 0; j < i + 1; j++) {
					pontosDif[j] = pontos.get(j);
				}

				double diferenca = difDivididas(pontosDif);

				for (int w = 0; w < i; w++) {
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
