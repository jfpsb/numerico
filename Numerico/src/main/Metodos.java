package main;

public class Metodos {

	public static void main(String[] args) {
		double a = 1;
		double b = 6;
		double precisao1 = 0.000001, precisao2 = precisao1;

		bissecao(a, b, precisao1);

		posicaoFalsa(a, b, precisao1, precisao2);

		pontoFixo(1.5, precisao1, precisao2);
		
		newtonRaphson(1.5, precisao1, precisao2);
		
		secante(1.5, 1.7, precisao1, precisao2);
	}

	private static void bissecao(double a, double b, double precisao) {
		double iteracoes = 0;
		double x = 0;

		if (possuiRaiz(a, b)) {
			if ((b - a) < precisao) {
				// Se intervalo for menor que precisão qualquer valor dentro do intervalo é raiz
				System.out.printf("\\nBISSEÇÃO >>> Qualquer valor no intervalo [%.3f, %.3f] com precisao %.3f é raiz!", a, b, precisao);
			} else {
				// Calcula número de iterações
				iteracoes = Math.ceil((Math.log10(b - a) - Math.log10(precisao)) / Math.log10(2));

				for (int i = 0; i < iteracoes; i++) {
					// Calcula ponto médio
					x = (a + b) / 2;

					if (funcao(x) == 0) {
						// Se ponto médio for raiz loop termina;
						break;
					}

					if ((funcao(a) * funcao(x)) < 0) {
						// Raiz está na metade à esquerda
						b = x;
					} else {
						// Raiz está na metade à direita
						a = x;
					}
				}

				System.out.printf("\nBISSEÇÃO >>> Raiz aproximada é x = %.3f", x);
			}
		} else {
			System.out.println("Não há raiz no intervalo!");
		}
	}

	private static void posicaoFalsa(double a, double b, double precisao1, double precisao2) {
		double iteracoes = 0;
		double x = 0;

		if (possuiRaiz(a, b)) {
			if ((b - a) < precisao1) {
				// Se intervalo for menor que precisão qualquer valor dentro do intervalo é raiz
				System.out.printf("\\nPOSIÇÃO FALSA >>> Qualquer valor no intervalo [%.3f, %.3f] com precisao %.3f é raiz!", a, b,
						precisao1);
			} else if (Math.abs(funcao(a)) < precisao2 || Math.abs(funcao(b)) < precisao2) {
				System.out.printf("A raiz pode ser %.3f ou %.3f", a, b);
			} else {
				// Calcula número de iterações
				iteracoes = Math.ceil((Math.log10(b - a) - Math.log10(precisao1)) / Math.log10(2));

				for (int i = 0; i < iteracoes; i++) {
					x = ((a * funcao(b)) - (b * funcao(a))) / (funcao(b) - funcao(a));

					if (funcao(x) == 0) {
						// Se ponto médio for raiz loop termina;
						break;
					}

					if (Math.abs(funcao(x)) < precisao2) {
						break;
					}

					if ((funcao(a) * funcao(x)) < 0) {
						// Raiz está na metade à esquerda
						b = x;
					} else {
						// Raiz está na metade à direita
						a = x;
					}
				}

				System.out.printf("\nPOSIÇÃO FALSA >>> Raiz aproximada é x = %.3f", x);
			}
		} else {
			System.out.println("Não há raiz no intervalo!");
		}
	}

	private static void pontoFixo(double aproximacao, double precisao1, double precisao2) {
		double x = 0;

		if (Math.abs(funcao(aproximacao)) < precisao1) {
			System.out.printf("PONTO FIXO >>> A raiz aproximada é igual a aproximação: %.3f", aproximacao);
		} else {
			while (true) {
				x = funcaoIteracao(aproximacao);

				if (Math.abs(funcao(x)) < precisao1 || Math.abs(x - aproximacao) < precisao2) {
					break;
				}

				aproximacao = x;
			}

			System.out.printf("\nPONTO FIXO >>> A raiz aproximada é igual a: %.5f", x);
		}
	}

	private static void newtonRaphson(double aproximacao, double precisao1, double precisao2) {
		double x = 0;

		if (Math.abs(funcao(aproximacao)) < precisao1) {
			System.out.printf("\\nNEWTON-RAPHSON >>> A raiz aproximada é igual a aproximação: %.3f", aproximacao);
		} else {
			while (true) {
				x = aproximacao - (funcao(aproximacao)/funcaoDerivada(aproximacao));

				if (Math.abs(funcao(x)) < precisao1 || Math.abs(x - aproximacao) < precisao2) {
					break;
				}

				aproximacao = x;
			}

			System.out.printf("\nNEWTON-RAPHSON >>> A raiz aproximada é igual a: %.5f", x);
		}
	}
	
	private static void secante(double aproximacao1, double aproximacao2, double precisao1, double precisao2) {
		double x = 0;
		if (Math.abs(funcao(aproximacao1)) < precisao1) {
			System.out.printf("\\nSECANTE >>> A raiz aproximada é igual a aproximação 1: %.3f", aproximacao1);
		} else if(Math.abs(funcao(aproximacao2)) < precisao1 || Math.abs(aproximacao2 - aproximacao1) < precisao2){
			System.out.printf("\\nSECANTE >>> A raiz aproximada é igual a aproximação 2: %.3f", aproximacao2);
		} else {
			while(true) {
				x = aproximacao2 - (funcao(aproximacao2)/(funcao(aproximacao2) - funcao(aproximacao1)))*(aproximacao2 - aproximacao1);
				
				if(Math.abs(funcao(x)) < precisao1 || Math.abs(x - aproximacao2) < precisao2) {
					break;
				}
				
				aproximacao1 = aproximacao2;
				aproximacao2 = x;
			}
			
			System.out.printf("\nSECANTE >>> A raiz aproximada é igual a: %.5f", x);
		}
	}

	private static double funcao(double x) {
		return Math.pow(x, 2) + x - 6;
	}

	private static double funcaoIteracao(double x) {
		return Math.sqrt(6 - x);
	}

	private static double funcaoDerivada(double x) {
		return (2 * x) + 1;
	}

	private static boolean possuiRaiz(double a, double b) {
		if ((funcao(a) * funcao(b)) < 0)
			return true;

		return false;
	}

}
