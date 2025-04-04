"""
Lumache - Python library for cooks and food lovers.
"""
from typing import Optional

__version__ = "0.1.0"


class InvalidKindError(Exception):
    """
    Raised if the kind is invalid.
        testing indentation.
    **tesitng**  
    for line change \n

    how it should work 
    or not  
    or how to elegantly
    
    gene_csv:
        * Assuming rows are genes, columns are samples (or cells)
            *testing
    """
    pass


def get_random_ingredients(kind:Optional[list[str]]=None):
    """
    Return a list of random ingredients as strings.
        *testing
    
    :param kind: (Optional[list[str]]) Optional "kind" of ingredients.
    :raise lumache.InvalidKindError: If the kind is invalid.
    :return: The ingredients list.
    :rtype: list[str]
    """
    return ["shells", "gorgonzola", "parsley"]
