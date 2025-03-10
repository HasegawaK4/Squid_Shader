using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;

public class GaussianBlur : CustomPostProcessVolumeComponent, IPostProcessComponent
{
    public ClampedFloatParameter blurSize = new ClampedFloatParameter(1f, 0f, 10f);

    Material material;

    public override void Setup()
    {
        if (material == null)
        {
            Shader shader = Shader.Find("Hidden/Custom/BlurPostProcess");
            material = new Material(shader);
        }
    }

    public override void Render(CommandBuffer cmd, HDCamera camera, RTHandle source, RTHandle destination)
    {
        material.SetFloat("_BlurSize", blurSize.value);
        HDUtils.DrawFullScreen(cmd, material, destination, source);
    }

    public bool IsActive() => blurSize.value > 0f;
}