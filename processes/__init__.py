from pathlib import Path


def model_path_resolution(model_source: str) -> str:
    """Resolve a model path relative to the project root."""
    if not model_source.startswith(('http://', 'https://')):
        model_path = Path(model_source)
        if not model_path.is_absolute():
            project_root = Path(__file__).parent.parent
            model_path = project_root / model_path
        model_source = str(model_path)
    return model_source
