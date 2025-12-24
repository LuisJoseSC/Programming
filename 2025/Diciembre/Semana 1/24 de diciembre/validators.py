import math
from typing import Tuple

def validate_coefficients(a,b,c) -> Tuple [float, float, float]:
    try:
        a = float(a)
        b = float(b)
        c = float(c)
    except (TypeError, ValueError):
        raise ValueError("Cofficients must be numeric values.")
    
    if not all(math.isfinite(x) for x in (a,b,c)):
        raise ValueError("Coffecients must be finite numbers.")
    if a == 0:
        raise ValueError("Coefficien 'a' must be non-zero for a quadratic equation. ")
    return a,b,c