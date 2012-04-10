from distutils.core import setup
setup(name="OASIS",
      author="David Robinson",
      author_email="dgrtwo@princeton.edu",
      description="Optimized Annotation System for Insertion Sequences",
      version="1.0",
      packages=["OASIS"],
      scripts=["scripts/OASIS"]
      requires=["biopython (>=1.46)"],
      )
