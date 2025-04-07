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
            
        <detach_from> options:
            - 'pre_encoder' (lat_pre will be detached, pre_encoder will not be trained);
            - 'encoder' (main_lat, pheno_lat, signature_lat will be detached, neither pre-encoder nor encoder will be trained)
        
        <pheno_meta>: For more details of the JSON structure, see :func:`utils.data_transformations`.        
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
