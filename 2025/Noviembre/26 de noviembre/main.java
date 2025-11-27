import java.util.Scanner;

public class main {

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        System.out.println("Calculadora de mierda pq java es una mierda.");
        System.out.println("1. Suma");
        System.out.println("2. Resta");
        System.out.println("3. Multiplicación");
        System.out.println("4. División");
        System.out.print("Seleccione la operación: ");
        int op = sc.nextInt();

        System.out.print("Ingrese el primer número: ");
        int a = sc.nextInt();

        System.out.print("Ingrese el segundo número: ");
        int b = sc.nextInt();

        int resultado = 0;

        switch (op) {
            case 1:
                resultado = suma.sumar(a, b);
                break;
            case 2:
                resultado = resta.restar(a, b);
                break;
            case 3:
                resultado = multiplicacion.multiplicar(a, b);
                break;
            case 4:
                resultado = division.dividir(a, b);
                break;
            default:
                System.out.println("Opción inválida.");
                return;
        }

        System.out.println("Resultado: " + resultado);
    }
}
