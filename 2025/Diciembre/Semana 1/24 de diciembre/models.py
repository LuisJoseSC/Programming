from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass
class QuadratricResult:
    a: float
    b: float
    c: float
    discriminant: float
    root_type: str
    roots: Optional[Tuple[float, ...]]
    factorization: Optional[str]