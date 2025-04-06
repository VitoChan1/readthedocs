"""
Lumache - Python library for cooks and food lovers.
"""
from typing import Optional

__version__ = "0.1.0"


class InvalidKindError(Exception):
    """
    Raised if the kind is invalid.
        testing indentation.
    .. note::
        For more details of the transformations, see :func:`utils.data_transformations`.

        <gene_meta> example:
        .. code-block::
        
            {
                'all': {
                    'gene_list': "*",
                    'pre_procedure': [],
                    'post_procedure': [{
                    'type': "ToTensor"
                    }]
                }
            }
            
        <pheno_meta>: For more details of the JSON structure, see :func:`utils.data_transformations`.

        <na_filter>: For phenotype data without any NA values, passing <na_filter>=False can
        improve the performance of reading a large file.
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
