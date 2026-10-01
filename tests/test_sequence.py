from app.sequence import parse_sequence, gc_percent, revcomp

def test_parse():
    h,s=parse_sequence('>x\nATGC U')
    assert h=='x' and s=='ATGCT'

def test_gc():
    assert gc_percent('ATGC')==50.0

def test_revcomp():
    assert revcomp('ATGC')=='GCAT'
