#Autor: Luis Jose Sanchez Carreño
#Fecha: Octubre 23 de 2025
#Proyecto diario.

#Calculadora de operaciones sencillas.

#Importaciones
import Division as D
import Multiplicacion as M
import Resta as R
import Suma as S

#Logica central
while True:
    print("Calculadora sencilla")
    print("Autor: Luis Jose Sanchez Carreño")
    print("1. Suma.")
    print("2. Resta.")
    print("3. Multiplicacion.")
    print("4. Division.")
    print("5. Modulo.")
    print("6. Potencia.")
    print("0. Salir.")
    opcion = int(input("Ingrese la operacion que quiere usar: "))
    
    if opcion == 1:
        num_1 = int(input("Ingrese el numero 1: "))
        num_2 = int(input("Ingrese el numero 2: "))
        resultado = S.suma(num_1, num_2)
        print(resultado)
    elif opcion == 2:
        num_1 = int(input("Ingrese el numero 1: "))
        num_2 = int(input("Ingrese el numero 2: "))
        resultado = R.resta(num_1, num_2)
        print(resultado)
    elif opcion == 3:
        num_1 = int(input("Ingrese el numero 1: "))
        num_2 = int(input("Ingrese el numero 2: "))
        resultado = M.multiplicacion(num_1, num_2)
        print(resultado)
    elif opcion == 4:
        num_1 = int(input("Ingrese el numero 1: "))
        num_2 = int(input("Ingrese el numero 2: "))
        resultado = D.division(num_1, num_2)
        print(resultado)
    elif opcion == 5:
        num_1 = int(input("Ingrese el numero 1: "))
        num_2 = int(input("Ingrese el numero 2: "))
        resultado = S.suma(num_1, num_2)
        print(resultado)
    elif opcion == 6:
        num_1 = int(input("Ingrese el numero 1: "))
        num_2 = int(input("Ingrese el numero 2: "))
        resultado = S.suma(num_1, num_2)
        print(resultado)
    elif opcion == 0:
        print("Gracias por usar la aplicacion.")
        break