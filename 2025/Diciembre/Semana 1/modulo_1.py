#Discriminante de una eucacion cuadratica.

def discriminante(a: int, b: int, c: int)-> str:
    if a == 0:
        raise ValueError("No es una ecuacion cuadratica")
    discriminante = b**2-4*a*c
    
    if discriminante > 0:
        return "tiene dos raices"
    elif discriminante == 0:
        return "tiene una raiz"
    else:
        return "sin raices reales"

print(discriminante(1,4,5))