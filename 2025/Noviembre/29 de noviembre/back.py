#Back end del proyecto.

empleados = [
    {"id": 1, "nombre": "Ana", "salario": 3200000},
    {"id": 2, "nombre": "Luis", "salario": 5400000},
    {"id": 3, "nombre": "Carlos", "salario": 2800000},
    {"id": 4, "nombre": "María", "salario": 4100000},
    {"id": 5, "nombre": "Sofía", "salario": 3600000}
]

def ordenar_salario(empleados):
    empleados_ordenados = sorted(empleados, key=lambda emp: emp["salario"], reverse=True)
    return empleados_ordenados

resultado = ordenar_salario(empleados)
print(resultado)
