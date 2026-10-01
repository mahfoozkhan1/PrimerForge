from typing import Optional, Literal, Any
from pydantic import BaseModel, Field

AssayType = Literal["pcr", "rtpcr", "qpcr", "ddpcr", "mirna", "sirna"]

class AnalyzeRequest(BaseModel):
    assay_type: AssayType
    sequence: str = Field(min_length=1)
    organism: Optional[str] = None
    run_blast: bool = False
    blast_db: str = "nt"
    max_pairs: int = Field(default=10, ge=1, le=50)
    custom_params: dict[str, Any] = Field(default_factory=dict)

class AnalyzeResponse(BaseModel):
    status: str
    assay_type: str
    sequence_length: int
    candidates: list[dict[str, Any]]
    best: Optional[dict[str, Any]] = None
    warnings: list[str] = Field(default_factory=list)
